using LinearSolve 
using SparseArrays
using Printf
using CairoMakie


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

taud_acx = fill(0.0,nx,ny);
taud_acy = fill(0.0,nx,ny);

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
        
        H_mid = 0.5*(H[i,j]+H[ip1,j]);
        taud_acx[i,j] = an["ρ"] * an["g"] * H_mid * (z_srf[ip1,j]-z_srf[i,j])/dx;

        H_mid = 0.5*(H[i,j]+H[i,jp1]);
        taud_acy[i,j] = an["ρ"] * an["g"] * H_mid * (z_srf[i,jp1]-z_srf[i,j])/dx;

    end
end

# Hack to ensure driving stress is ok in last grid point
taud_acx[end,:] = taud_acx[end-1,:];

# Function to calculate velocity 

function ij2n_ux(i,j,nx,ny)

    n = (i-1)*ny + j

    return n 
end

function ij2n_uy(i,j,nx,ny)

    n = (i-1)*ny + j + nx*ny
    
    return n 
end

function calc_vel_ssa!(ux,uy,H,μ,taud_acx,taud_acy,β_acx,β_acy,dx)
    # Calculate the diagnostic SSA velocity solution 
    # given ice thickness, viscosity, driving stress and basal friction coefficient

    # Get some constants 
    nx, ny = size(H);

    dy = dx;

    inv_dxdx    = 1.0 / (dx*dx);
    inv_dydy    = 1.0 / (dy*dy);
    inv_dxdy    = 1.0 / (dx*dy);

    # Define vertically-integrated viscosity 
    N_aa = H .* μ;

    # Stagger N to ab-nodes
    N_ab = copy(N_aa);

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
            
            N_ab[i,j] = 0.25 * (N_aa[i,j]+N_aa[ip1,j]+N_aa[i,jp1]+N_aa[ip1,jp1]);

        end 
    end 

    #
    #
    # Populate SSA stress balance matrix equation Ax = b 
    # [ A_ux   A_vx ]  [ u ]  = [ b_x ]
    # [ A_uy   A_vy ]  [ v ]    [ b_y ]
    #
    #

    n_terms = 9;
    n_u     = 2*nx*ny;
    n_sprs  = n_u*n_terms;

    # Populate dense array for now
    u = fill(0.0,n_u);
    b = fill(0.0,n_u);

    use_dense = true;

    if use_dense
        A = fill(0.0,n_u,n_u);
    else
        Ai = fill(0,  n_sprs);
        Aj = fill(0,  n_sprs);
        Av = fill(0.0,n_sprs);
    end

    # Equation is being defined for acx-nodes (x-direction equation)

    k = 0 

    for i in 1:nx 
        for j in 1:ny 

            # Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            # Set the row in matrix A that the equation is being defined for:
            nr = (i-1)*ny + j

if use_dense
            # -- vx terms -- 
            
            # ux(i,j)
            n = ij2n_ux(i,j,nx,ny);
            A[nr,n] =   -4.0*inv_dxdx*N_aa[ip1,j]
                        -4.0*inv_dxdx*N_aa[i,j]
                        -1.0*inv_dydy*N_ab[i,j]
                        -1.0*inv_dydy*N_ab[i,jm1]
                        -β_acx[i,j];
            
            # ux(i+1,j)
            n = ij2n_ux(ip1,j,nx,ny);
            A[nr,n] =    4.0*inv_dxdx*N_aa[ip1,j];
            
            # ux(i-1,j)
            n = ij2n_ux(im1,j,nx,ny);
            A[nr,n] =    4.0*inv_dxdx*N_aa[i,j];

            # ux(i,j+1)
            n = ij2n_ux(i,jp1,nx,ny);
            A[nr,n] =    1.0*inv_dydy*N_ab[i,j];

            # ux(i,j-1)
            n = ij2n_ux(i,jm1,nx,ny);
            A[nr,n] =    1.0*inv_dydy*N_ab[i,jm1];

            # -- vy terms -- 

            # uy(i,j)
            n = ij2n_uy(i,j,nx,ny);
            A[nr,n] =   -2.0*inv_dxdy*N_aa[i,j]
                        -1.0*inv_dxdy*N_ab[i,j];
            
            # uy(i+1,j)
            n = ij2n_uy(ip1,j,nx,ny);
            A[nr,n] =    2.0*inv_dxdy*N_aa[ip1,j]
                        +1.0*inv_dxdy*N_ab[i,j];
            
            # uy(i+1,j-1)
            n = ij2n_uy(ip1,jm1,nx,ny);
            A[nr,n] =   -2.0*inv_dxdy*N_aa[ip1,j]
                        -1.0*inv_dxdy*N_ab[i,jm1];
            
            # uy(i,j-1)
            n = ij2n_uy(i,jm1,nx,ny);
            A[nr,n] =    2.0*inv_dxdy*N_aa[i,j]
                        +1.0*inv_dxdy*N_ab[i,jm1];

else
            # -- vx terms -- 
            
            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = -4.0*inv_dxdx*N_aa[ip1,j]
                    -4.0*inv_dxdx*N_aa[i,j]
                    -1.0*inv_dydy*N_ab[i,j]
                    -1.0*inv_dydy*N_ab[i,jm1]
                    -β_acx[i,j];

            # ux(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(ip1,j,nx,ny);
            Av[k] =  4.0*inv_dxdx*N_aa[ip1,j];

            # ux(i-1,j)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] =  4.0*inv_dxdx*N_aa[i,j];

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] =  1.0*inv_dydy*N_ab[i,j];

            # ux(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jm1,nx,ny);
            Av[k] =  1.0*inv_dydy*N_ab[i,jm1];
            
            # -- vy terms -- 

            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = -2.0*inv_dxdy*N_aa[i,j]
                    -1.0*inv_dxdy*N_ab[i,j];
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = -2.0*inv_dxdy*N_aa[ip1,j]
                    +1.0*inv_dxdy*N_ab[i,j];
            
            # uy(i+1,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,jm1,nx,ny);
            Av[k] = -2.0*inv_dxdy*N_aa[ip1,j]
                    -1.0*inv_dxdy*N_ab[i,jm1];
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] =  2.0*inv_dxdy*N_aa[i,j]
                    +1.0*inv_dxdy*N_ab[i,jm1];
             
end

            # [u] value
            u[nr] = ux[i,j];

            # [b] value 
            b[nr] = taud_acx[i,j];

        end
    end

    # Equation is being defined for acy-nodes (y-direction equation)

    for i in 1:nx 
        for j in 1:ny 

            # Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            # Set the row in matrix A that the equation is being defined for:
            nr = (i-1)*ny + j + nx*ny

if use_dense
            # -- vy terms -- 
            
            # uy(i,j)
            n = ij2n_uy(i,j,nx,ny);
            A[nr,n] =   -4.0*inv_dydy*N_aa[i,jp1]
                        -4.0*inv_dydy*N_aa[i,j]
                        -1.0*inv_dxdx*N_ab[i,j]
                        -1.0*inv_dxdx*N_ab[im1,j]
                        -β_acy[i,j];

            # uy(i,j+1)
            n = ij2n_uy(i,jp1,nx,ny);
            A[nr,n] =    4.0*inv_dydy*N_aa[i,jp1];
            
            # uy(i,j-1)
            n = ij2n_uy(i,jm1,nx,ny);
            A[nr,n] =    4.0*inv_dydy*N_aa[i,j];
            
            # uy(i+1,j)
            n = ij2n_uy(ip1,j,nx,ny);
            A[nr,n] =    1.0*inv_dxdx*N_ab[i,j];
            
            # uy(i-1,j)
            n = ij2n_uy(im1,j,nx,ny);
            A[nr,n] =    1.0*inv_dxdx*N_ab[im1,j];
            
            # -- vx terms -- 

            # ux(i,j)
            n = ij2n_ux(i,j,nx,ny);
            A[nr,n] =   -2.0*inv_dxdy*N_aa[i,j]
                        -1.0*inv_dxdy*N_ab[i,j];

            # ux(i,j+1)
            n = ij2n_ux(i,jp1,nx,ny);
            A[nr,n] =    2.0*inv_dxdy*N_aa[i,jp1]
                        +1.0*inv_dxdy*N_ab[i,j];

            # ux(i-1,j+1)
            n = ij2n_ux(im1,jp1,nx,ny);
            A[nr,n] =   -2.0*inv_dxdy*N_aa[i,jp1]
                        -1.0*inv_dxdy*N_ab[im1,j];

            # ux(i-1,j)
            n = ij2n_ux(im1,j,nx,ny);
            A[nr,n] =    2.0*inv_dxdy*N_aa[i,j]
                        +1.0*inv_dxdy*N_ab[im1,j];
else

            # -- vy terms -- 
            
            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = -4.0*inv_dydy*N_aa[i,jp1]
                    -4.0*inv_dydy*N_aa[i,j]
                    -1.0*inv_dxdx*N_ab[i,j]
                    -1.0*inv_dxdx*N_ab[im1,j]
                    -β_acy[i,j]; 

            # uy(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jp1,nx,ny);
            Av[k] =  4.0*inv_dydy*N_aa[i,jp1];  
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] =  4.0*inv_dydy*N_aa[i,j];    
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] =  1.0*inv_dxdx*N_ab[i,j];     
            
            # uy(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(im1,j,nx,ny);
            Av[k] =  1.0*inv_dxdx*N_ab[im1,j];   
            
            # -- vx terms -- 

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = -2.0*inv_dxdy*N_aa[i,j]
                    -1.0*inv_dxdy*N_ab[i,j]; 

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] =  2.0*inv_dxdy*N_aa[i,jp1]
                    +1.0*inv_dxdy*N_ab[i,j];

            # ux(i-1,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,jp1,nx,ny);
            Av[k] = -2.0*inv_dxdy*N_aa[i,jp1]
                    -1.0*inv_dxdy*N_ab[im1,j];
            
            # ux(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] =  2.0*inv_dxdy*N_aa[i,j]
                    +1.0*inv_dxdy*N_ab[im1,j];
            
end
            # [u] value
            u[nr] = uy[i,j];

            # [b] value 
            b[nr] = taud_acy[i,j];
            
        end
    end

    # Now A (dense), x and b have been populated

    if use_dense
        Asp = sparse(A);
    else
        Asp = sparse(Ai,Aj,Av);
    end

    use_linsolve = true

    if use_linsolve
        prob = LinearProblem(Asp, b; u0=u);
        sol = solve(prob);
        unew = sol.u;

    else

        unew = Asp \ b;

    end

    for i = 1:nx
        for j = 1:ny
            n = ij2n_ux(i,j,nx,ny);
            ux[i,j] = unew[n];
        end
    end
    for i = 1:nx
        for j = 1:ny
            n = ij2n_uy(i,j,nx,ny);
            uy[i,j] = unew[n];
        end
    end

    return Asp, u, b
end

function plot_out(var)

    fig,ax,hm = heatmap(var)
    Colorbar(fig[1,end+1],hm)
    save("test.pdf",fig)

    println("extrema: ",extrema(var))
end


# Test 
A, u, b = calc_vel_ssa!(ux,uy,H,μ,taud_acx,taud_acy,β_acx,β_acy,dx);
println("ux, uy: ", extrema(ux), " | ", extrema(uy) )

plot_out(ux)