using NCDatasets
using Printf
using CairoMakie

include(pwd()*"/src/pagos-base.jl");

"""
    calc_vel_stream!(p)

    Calculate analytical solution of stream 
"""
function define_stream_slab(;H0=1000.0,μ0=1e5,β0=1e3,α=1e-2,ρ=910.0,g=9.81,verbose=true)

    # Get constants
    η  = β0 * H0 / μ0;
    F2 = H0 / (3.0 * μ0);

    # Calculate rate factor
    # ATT = (2.0*visc_eff)^(-1) 

    # Calculate analytical velocity solution
    ub = ρ * g * H0 * α / β0; 
    u  = ub + ρ * g * H0^2 * α / (3.0 * μ0);

    # Define dictionary with values of analytical stream solution
    strm = Dict("H0"=>H0,"μ0"=>μ0,"β0"=>β0,"α"=>α,"ρ"=>ρ,"η"=>η,"g"=>g,"F2"=>F2,"ub"=>ub,"u"=>u)
    
    # Print summary
    if verbose
        println("Analytical stream defined:")
        for (key,val) in strm
            @printf("  %4s = %6.3g",key,val)
        end
        print("\n")
    end

    return strm
end

@inline function schoof2006_vel(y,taud,H0,B,L,W,m)

    # Calculate the analytical solution for u
    ua = -2.0 * taud^3 * L^4 / (B^3 * H0^3)
    ub = ( 1.0 / 4.0                      ) * (   (y/L)^     4.0  - (m+1.0)^(     4.0/m ))
    uc = (-3.0 / ((m+1.0)   * (    m+4.0))) * (abs(y/L)^(  m+4.0) - (m+1.0)^(1.0+(4.0/m)))
    ud = ( 3.0 / ((m+1.0)^2 * (2.0*m+4.0))) * (abs(y/L)^(2*m+4.0) - (m+1.0)^(2.0+(4.0/m)))
    ue = (-1.0 / ((m+1.0)^3 * (3.0*m+4.0))) * (abs(y/L)^(3*m+4.0) - (m+1.0)^(3.0+(4.0/m)))
    u = ua * (ub + uc + ud + ue);
    
    # Outside the ice-stream, velocity is zero
    if abs(y) > W || abs(u) < 1e-10
        u = 0.0
    end

    return u
end

function define_stream_schoof2006(dx;xmax=140e3,ymax=50e3,H0=1e3,rf=1e-16,α=1e-3,ρ=910.0,g=9.81,n_glen=3,W=25e3,m=1.55)

    # Schoof (2006) domain - constant slope slab

    # Intialize domain 

    # Define axes x,y ###

    xmin = 0.0;
    ymin = -ymax;

    nx = Int(xmax/dx)+1;
    ny = Int((ymax-ymin)/dx)+1;
    
    x0 = 0.0;
    y0 = -(ny-1)/2*dx;

    xc = [xmin + (i-1)*dx for i in 1:nx];
    yc = [ymin + (i-1)*dx for i in 1:ny];

    # ===== Intialize topography and set parameters =========

    H = fill(H0,nx,ny);
    
    z_bed = fill(NaN,nx,ny);

    for i = 1:nx
        z_bed[i,:] .= 10000.0 .- α .* (xc[i] .- xmin);
    end

    # Define surface elevation 
    z_srf = z_bed .+ H;

    # Calculate the gravitational driving stress f
    taud = ρ * g * H0 * α;
    
    # Calculate the ice hardness factor B
    B = rf^(-1/n_glen);
    
    # Determine constant L (ice-stream width)
    L = W / ((1.0+m)^(1.0/m));

    # Calculate the till yield stress across the stream
    # and analytical velocity solution
    tau_c = fill(NaN,size(H));
    ux    = fill(NaN,nx,ny);
    
    for j = eachindex(yc)
        tau_c[:,j] .= taud * abs(yc[j] / L)^m;
        ux[:,j]    .= schoof2006_vel(yc[j],taud,H0,B,L,W,m);
    end

    
    # Assign analytical values (tau_c as a boundary condition, ux as initial condition)
    cb = tau_c;

    # Determine constant L too, for diagnostic output
    L = W / ((1.0+m)^(1.0/m));

    println("SLAB-S06: H0      = ", H0)
    println("SLAB-S06: alpha   = ", α)
    println("SLAB-S06: W       = ", W)
    println("SLAB-S06: L       = ", L) 
    println("SLAB-S06: m       = ", m) 
    println("SLAB-S06: ρ g     = ", ρ, " ", g)
    println("SLAB-S06: f       = ", (ρ*g*H0)*α)
    println("SLAB-S06: ATT     = ", rf)
    println("SLAB-S06: cb      = ", cb[1,1])
    println("SLAB-S06: tau_c   = ", tau_c[1,1])
    println("SLAB-S06: max(ux) = ", maximum(ux))

    p = Dict("H0"=>H0,"α"=>α,"W"=>W,"L"=>L,"m"=>m,"ρ"=>ρ,"g"=>g,"rf"=>rf);

    return Dict("p"=>p, "xc"=>xc, "yc"=>yc, "H"=>H, "z_bed"=>z_bed, "z_srf"=>z_srf, "ux"=>ux, "tau_c"=>tau_c, "cb"=>cb)
end

function solve_stream_slab(an)

    # Define axes x,y ###

    xmin = 0.0;
    ymin = -(ny-1)/2*dx;

    xc = [xmin + (i-1)*dx for i in 1:nx];
    yc = [ymin + (i-1)*dx for i in 1:ny];

    # Define variables with initial values ###

    # Ice thickness
    H = fill(an["H0"],nx,ny); 

    # Bed elevation
    z_bed = fill(0.0,nx,ny);
    for j in 1:ny
        z_bed[:,j] = [10000.0 - an["α"]*(x) for x in xc];
    end

    # Surface elevation
    z_srf = z_bed .+ H;

    # Viscosity
    μ = fill(an["μ0"],nx,ny);

    # Basal friction 
    β = fill(an["β0"],nx,ny);
    β_acx, β_acy = stagger_beta(β);

    # Get driving stress
    taud_acx, taud_acy = calc_driving_stress(H,z_srf,dx,dx,an["ρ"],an["g"]);

    # Now, solve for velocity:
    ux = fill(0.0,nx,ny);
    uy = fill(0.0,nx,ny);
    calc_vel_ssa!(ux,uy,H,μ,taud_acx,taud_acy,β_acx,β_acy,dx);
    
    #println("ux, uy: ", extrema(ux), " | ", extrema(uy) )

    # Collect variables for output 

    return Dict("xc"=>xc,"yc"=>yc,"H"=>H,"z_bed"=>z_bed,"z_srf"=>z_srf,"μ"=>μ,"β"=>β,"β_acx"=>β_acx,"β_acy"=>β_acy,
                    "taud_acx"=>taud_acx,"taud_acy"=>taud_acy,"ux"=>ux,"uy"=>uy)
end

function plot_out(var)

    fig,ax,hm = heatmap(var)
    Colorbar(fig[1,end+1],hm)
    save("test.pdf",fig)

    println("extrema: ",extrema(var))
end


######################################
# Parameters #

ρ  = 910.0
g  = 9.81

nx = 11
ny = 3
dx = 5e3

# Analytical cases
an1 = define_stream_slab(H0=1000.0,μ0=1e5,β0=1e3, α=1e-3,ρ=ρ,g=g);
an2 = define_stream_slab(H0= 500.0,μ0=4e5,β0=30.0,α=1e-3,ρ=ρ,g=g);

######################################

# Test
strm1 = solve_stream_slab(an1);
strm2 = solve_stream_slab(an2);

plot_out(strm1["ux"])
plot_out(strm2["ux"])







########################################################################

# To do: 1D numerical case translated from Matlab:

function calc_vel_diva_1D!(ux,H,H0,mu,beta_sl,dx,ρ,g,α)
    #calc_vel_ssa!(ux,uy,H,μ,taud_acx,taud_acy,β_acx,dx)
    
    eta  = beta_sl*H0/mu; 
    
    rhs = ρ .* g .* (H[1:end-1] .+ H[2:end])./2 .* (diff(H)./dx .- α);
    F2  = (H[1:end-1] .+ H[2:end]) ./ 2 ./ (3*eta);
    B   = beta_sl ./ (1 .+ beta_sl .* F2);
   
    d0 = -B .- 4*eta*H[1:end-1]/dx^2 .- 4*eta*H[2:end]/dx^2;
    dr = 4 * eta * H[2:end] / dx^2;
    dl = 4 * eta * H[1:end-1] / dx^2;
    #A = spdiags([dr d0 dl],[-1 0 1],N,N)'; 
    A = spdiagm(-1 => dr, 0 => d0, 1 => dl);
    #A[1,N] = dl[1];
    #A[N,1] = dr[end];
   
    u = A\rhs;

    return
end


#ux1D = ux[:,2];
#H1D  = H[:,2];
#calc_vel_diva_1D!(ux1D[1:end-1],H1D,an["H0"],an["μ0"],an["β0"],dx,an["ρ"],an["g"],an["α"])
#println("ux1D: ", extrema(ux1D[1:end-1]))

println("Done.")