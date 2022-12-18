using NCDatasets
using Printf
using CairoMakie

include(pwd()*"/src/pagos-base.jl");

"""
    calc_vel_stream!(p)

    Calculate analytical solution of stream 
"""
function define_stream(;H0=1000.0,μ0=1e5,β0=1e3,α=1e-2,ρ=910.0,g=9.81,verbose=true)

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

    # Calculate analytical stream function to get tau_c and ux

    #ux, tau_c = SSA_Schoof2006_analytical_solution_yelmo(y,α,H0,rf_const,W,m,n_glen,ρ,g);
    ux    = fill(NaN,size(H));
    tau_c = fill(NaN,size(H));

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

    return p, xc, yc, H, z_bed, z_srf, ux, tau_c, cb
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
an1 = define_stream(H0=1000.0,μ0=1e5,β0=1e3, α=1e-3,ρ=ρ,g=g);
an2 = define_stream(H0= 500.0,μ0=4e5,β0=30.0,α=1e-3,ρ=ρ,g=g);

# Define analytical case of interest now 
an = an1; 

######################################

# Define axes x,y ###

x0 = 0.0;
y0 = -(ny-1)/2*dx;

xc = [x0 + (i-1)*dx for i in 1:nx];
yc = [y0 + (i-1)*dx for i in 1:ny];


# Define variables with initial values ###

# Ice thickness
H = fill(an["H0"],nx,ny); 

# Viscosity
μ = fill(an["μ0"],nx,ny);

# Basal friction 
β = fill(an["β0"],nx,ny);

# Stagger β to acx and acy nodes 
β_acx = copy(β);
β_acy = copy(β);

for i in 1:nx
    for j in 1:ny

        # Periodic boundary conditions in x and y 
        ip1 = i+1
        if ip1 == nx+1
            ip1 = 1 
        end
        jp1 = j+1
        if jp1 == ny+1
            jp1 = 1 
        end
        
        β_acx[i,j] = 0.5*(β[i,j]+β[ip1,j]);
        β_acy[i,j] = 0.5*(β[i,j]+β[i,jp1]);
    end
end

# Bed elevation
z_bed = fill(0.0,nx,ny);
for j in 1:ny
    z_bed[:,j] = [10000.0 - an["α"]*(x) for x in xc];
end

# Surface elevation
z_srf = z_bed .+ H;

@show z_srf[:,1]

ux = fill(0.0,nx,ny);
uy = fill(0.0,nx,ny);

taud_acx, taud_acy = calc_driving_stress(H,z_srf,dx,dx,an["ρ"],an["g"]);


# Now, solve for velocity:

A, u, b = calc_vel_ssa!(ux,uy,H,μ,taud_acx,taud_acy,β_acx,β_acy,dx);
#println("ux, uy: ", extrema(ux), " | ", extrema(uy) )

#ux1D = ux[:,2];
#H1D  = H[:,2];
#calc_vel_diva_1D!(ux1D[1:end-1],H1D,an["H0"],an["μ0"],an["β0"],dx,an["ρ"],an["g"],an["α"])
#println("ux1D: ", extrema(ux1D[1:end-1]))

plot_out(ux)


