using Parameters
using NCDatasets 

using LinearSolve 
using SparseArrays


mutable struct cell 
    x::Float32
    y::Float32 
    i::Int 
    j::Int
    ii::Array{Int}
    jj::Array{Int}
    mask::Int
end 

mutable struct LinearSolverClass

    nx::Integer
    ny::Integer
    n_terms::Integer
    nmax::Integer
    n_sprs::Integer
        
    n2i::Array{Int,1}
    n2j::Array{Int,1}
    ij2n::Array{Int,2}
        
    a_ptr::Array{Int,1}
    a_idx::Array{Int,1}

    A_val::Array{Real,1}
    b_val::Array{Real,1}
    x_val::Array{Real,1}

    LinearSolverClass(nx,ny,n_terms) = new(nx,ny,n_terms,
                                        2*nx*ny,
                                        2*nx*ny*n_terms,
                                        zeros(2*nx*ny),
                                        zeros(2*nx*ny),
                                        zeros(nx,ny),
                                        zeros(2*nx*ny*n_terms),
                                        zeros(2*nx*ny*n_terms),
                                        zeros(2*nx*ny*n_terms),
                                        zeros(nx*ny),
                                        zeros(nx,ny))

end

#mutable struct pagos

# Driving stress
function calc_driving_stress(H,z_srf,dx,dy,ρ,g)
    
    nx, ny = size(H);

    taud_acx = fill(0.0,nx,ny);
    taud_acy = fill(0.0,nx,ny);

    for i in 1:nx
        for j in 1:ny

            # BC: Periodic boundary conditions in x and y 
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1 
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1 
            end
            
            H_mid = 0.5*(H[i,j]+H[ip1,j]);
            taud_acx[i,j] = ρ * g * H_mid * (z_srf[ip1,j]-z_srf[i,j])/dx;

            H_mid = 0.5*(H[i,j]+H[i,jp1]);
            taud_acy[i,j] = ρ * g * H_mid * (z_srf[i,jp1]-z_srf[i,j])/dy;

        end
    end

    # BC: periodic...
    # Hack to ensure driving stress is ok in last grid point,
    # since z_srf might not be periodic...
    taud_acx[end,:] = taud_acx[end-1,:];

    return taud_acx, taud_acy
end


# Functions to calculate velocity 

function calc_F_integral(visc_eff,H_ice,f_ice,zeta_aa,n)

    # To do...
    
    return Fn
end

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

            # BC: Periodic boundary conditions in x and y 
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

    # Define array A component vectors (I, J, V)
    Ai = fill(0,  n_sprs);
    Aj = fill(0,  n_sprs);
    Av = fill(0.0,n_sprs);

    # Equation is being defined for acx-nodes (x-direction equation)

    k = 0 

    for i in 1:nx
        for j in 1:ny 

            # BC: Periodic boundary conditions
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
            
            # -- vx terms -- 
            
            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-4.0*inv_dxdx*N_aa[ip1,j]
                    -4.0*inv_dxdx*N_aa[i,j]
                    -1.0*inv_dydy*N_ab[i,j]
                    -1.0*inv_dydy*N_ab[i,jm1]
                    -β_acx[i,j]);

            # ux(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(ip1,j,nx,ny);
            Av[k] = ( 4.0*inv_dxdx*N_aa[ip1,j]);

            # ux(i-1,j)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 4.0*inv_dxdx*N_aa[i,j]);

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 1.0*inv_dydy*N_ab[i,j]);

            # ux(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jm1,nx,ny);
            Av[k] = ( 1.0*inv_dydy*N_ab[i,jm1]);
            
            # -- vy terms -- 

            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-2.0*inv_dxdy*N_aa[i,j]
                     -1.0*inv_dxdy*N_ab[i,j]);
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = (-2.0*inv_dxdy*N_aa[ip1,j]
                     +1.0*inv_dxdy*N_ab[i,j]);
            
            # uy(i+1,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,jm1,nx,ny);
            Av[k] = (-2.0*inv_dxdy*N_aa[ip1,j]
                     -1.0*inv_dxdy*N_ab[i,jm1]);
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 2.0*inv_dxdy*N_aa[i,j]
                     +1.0*inv_dxdy*N_ab[i,jm1]);

            # [u] value
            u[nr] = ux[i,j];

            # [b] value 
            b[nr] = taud_acx[i,j];

        end
    end

    # Equation is being defined for acy-nodes (y-direction equation)

    for i in 1:nx 
        for j in 1:ny

            # BC: Periodic boundary conditions
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

            # -- vy terms -- 
            
            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-4.0*inv_dydy*N_aa[i,jp1]
                     -4.0*inv_dydy*N_aa[i,j]
                     -1.0*inv_dxdx*N_ab[i,j]
                     -1.0*inv_dxdx*N_ab[im1,j]
                     -β_acy[i,j]); 

            # uy(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jp1,nx,ny);
            Av[k] = ( 4.0*inv_dydy*N_aa[i,jp1]);  
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 4.0*inv_dydy*N_aa[i,j]);    
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 1.0*inv_dxdx*N_ab[i,j]);     
            
            # uy(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(im1,j,nx,ny);
            Av[k] = ( 1.0*inv_dxdx*N_ab[im1,j]);   
            
            # -- vx terms -- 

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-2.0*inv_dxdy*N_aa[i,j]
                     -1.0*inv_dxdy*N_ab[i,j]); 

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 2.0*inv_dxdy*N_aa[i,jp1]
                     +1.0*inv_dxdy*N_ab[i,j]);

            # ux(i-1,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,jp1,nx,ny);
            Av[k] = (-2.0*inv_dxdy*N_aa[i,jp1]
                     -1.0*inv_dxdy*N_ab[im1,j]);
            
            # ux(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 2.0*inv_dxdy*N_aa[i,j]
                     +1.0*inv_dxdy*N_ab[im1,j]);
            
            # [u] value
            u[nr] = uy[i,j];

            # [b] value 
            b[nr] = taud_acy[i,j];
            
        end
    end

    # Now A (dense), x and b have been populated
    # Define sparse array A:

    # Now u, b and A components (I, J, V vectors) have been defined.
    # Convert into a sparse array for solving:

    Asp = sparse(Ai,Aj,Av);


    use_linsolve = true

    if use_linsolve
        prob = LinearProblem(Asp, b; u0=u);
        sol = solve(prob);
        unew = sol.u;

    else

        unew = Asp \ b;

    end

    # Update velocity arrays with new solution
    for i = 1:nx
        for j = 1:ny
            n = ij2n_ux(i,j,nx,ny);
            ux[i,j] = unew[n];
            n = ij2n_uy(i,j,nx,ny);
            uy[i,j] = unew[n];
        end
    end

    return Asp, u, b
end

println("Loaded pagos-base.jl.")


