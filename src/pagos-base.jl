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

function stagger_beta(β)

    nx, ny = size(β);

    # Stagger β to acx and acy nodes 
    β_acx = copy(β);
    β_acy = copy(β);

    for i in 1:nx
        for j in 1:ny

            # BC: Periodic boundary conditions
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

    return β_acx, β_acy
end

"""
    calc_beta_aa_power_plastic(ux_b,uy_b,c_bed,f_ice,q,u_0)

    Calculate basal friction coefficient (beta) that
    enters the SSA solver as a function of basal velocity
    using a power-law form following Bueler and van Pelt (2015)
    
"""
function calc_beta_aa_power_plastic(ux_b::Array{T,2},uy_b::Array{T,2},
                                    c_bed::Array{T,2},f_ice::Array{T,2},q,u_0) where T
    

    # Local variables
    ub_min    = 1e-3;               # [m/yr] Minimum velocity is positive small value to avoid divide by zero
    ub_sq_min = ub_min^2;

    nx, ny = size(ux_b);

    # Initially set friction to zero everywhere
    beta = fill(0.0,nx,ny);
    
    for i = 1:nx
        for j = 1:ny

            # Get neighbor indices 
            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1) 
            jp1 = min(j+1,ny)
            
            if f_ice[i,j] == 1.0 
                # Fully ice-covered point with some fully ice-covered neighbors 

                cb_aa = c_bed[i,j];

                if q == 1.0 
                    # Linear law, no f(ub) term

                    beta[i,j] = cb_aa / u_0;

                else
                    # Non-linear law with f(ub) term 

                    # Unstagger velocity components to aa-nodes 
                    ux_aa = 0.5*(ux_b[i,j]+ux_b[im1,j]);
                    uy_aa = 0.5*(uy_b[i,j]+uy_b[i,jm1]);
                    
                    uxy_aa = sqrt(ux_aa^2 + uy_aa^2 + ub_sq_min);
                    
                    if q == 0
                        # Plastic law

                        beta[i,j] = cb_aa * (1.0 / uxy_aa);
                        
                    else

                        beta[i,j] = cb_aa * (uxy_aa / u_0)^q * (1.0 / uxy_aa);
                    
                    end
                    
                end

            else
                # Assign minimum velocity value, no staggering for simplicity

                if q == 1.0 
                    # Linear law, no f(ub) term

                    beta[i,j] = c_bed[i,j] / u_0;

                else
                    
                    uxy_b  = ub_min;

                    if q == 0.0 
                        # Plastic law

                        beta[i,j] = c_bed[i,j] * (1.0 / uxy_b);

                    else

                        beta[i,j] = c_bed[i,j] * (uxy_b / u_0)^q * (1.0 / uxy_b);

                    end 

                end

            end

        end
    end

    return beta
    
end

function calc_visc_eff_2D_aa(ux,uy,ATT,H_ice,f_ice,dx,dy;n_glen=3,eps_0=1e-6)
    # Calculate 3D effective viscosity following L19, Eq. 2
    # Use of eps_0 ensures non-zero positive viscosity value everywhere 
    # Note: viscosity is first calculated on ab-nodes, then 
    # unstaggered back to aa-nodes. This ensures more stability for 
    # visc_eff (less likely to blow up for low strain rates). 

    visc_min = 1e5;

    nx, ny = size(ux);

    # Calculate exponents 
    p1 = (1.0 - n_glen)/(2.0*n_glen);
    p2 = -1.0/n_glen;

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0*eps_0;

    # Calculate visc_eff on aa-nodes

    visc = fill(visc_min,nx,ny);

    for i = 1:nx
        for j = 1:ny  

            if f_ice[i,j] == 1.0

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

                # Get strain rate terms
                dudx_aa = (ux[i,j]-ux[im1,j])/dx 
                dvdy_aa = (uy[i,j]-uy[i,jm1])/dy 
                
                dudy_aa_1 = (ux[i,jp1]-ux[i,jm1])/(2.0*dy)
                dudy_aa_2 = (ux[im1,jp1]-ux[im1,jm1])/(2.0*dy)
                dudy_aa   = 0.5*(dudy_aa_1+dudy_aa_2)

                dvdx_aa_1 = (uy[ip1,j]-uy[im1,j])/(2.0*dx)
                dvdx_aa_2 = (uy[ip1,jm1]-uy[im1,jm1])/(2.0*dx)
                dvdx_aa   = 0.5*(dvdx_aa_1+dvdx_aa_2)

                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_aa = dudx_aa^2 + dvdy_aa^2 + dudx_aa*dvdy_aa + 0.25*(dudy_aa+dvdx_aa)^2 + eps_0_sq

                # Get rate factor on central node
                ATT_aa = ATT[i,j];

                # Calculate effective viscosity on ab-nodes
                visc[i,j] = 0.5*(eps_sq_aa)^(p1) * ATT_aa^(p2)

            end

        end 
    end

    return visc

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

function calc_vel_ssa(ux,uy,H,μ,taud_acx,taud_acy,β_acx,β_acy,dx)
    # Calculate the diagnostic SSA velocity solution 
    # given ice thickness, viscosity, driving stress and basal friction coefficient

    use_sico_discretization = false 

    # Get some constants 
    nx, ny = size(H);

    dy = dx;
    
    dxdx = (dx*dx);
    dydy = (dy*dy);
    dxdy = (dx*dy);

    # Define vertically-integrated viscosity (aa-nodes)
    N = H .* μ;

    if use_sico_discretization

        # Stagger N to ab-nodes
        N_ab = copy(N);

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
                
                N_ab[i,j] = 0.25 * (N[i,j]+N[ip1,j]+N[i,jp1]+N[ip1,jp1]);

            end 
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

if use_sico_discretization
            
            # Standard (sicopolis) discretization
            
            # -- vx terms -- 
            
            # ux(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(ip1,j,nx,ny);
            Av[k] = ( 4.0/dxdx*N[ip1,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-4.0/dxdx * (N[ip1,j]+N[i,j])
                     -1.0/dydy * (N_ab[i,j]+N_ab[i,jm1])
                    -β_acx[i,j]);

            # ux(i-1,j)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 4.0/dxdx*N[i,j]);

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 1.0/dydy*N_ab[i,j]);

            # ux(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jm1,nx,ny);
            Av[k] = ( 1.0/dydy*N_ab[i,jm1]);
            
            # -- vy terms -- 

            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,j]
                     -1.0/dxdy*N_ab[i,j]);
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[ip1,j]
                     +1.0/dxdy*N_ab[i,j]);
            
            # uy(i+1,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,jm1,nx,ny);
            Av[k] = (-2.0/dxdy*N[ip1,j]
                     -1.0/dxdy*N_ab[i,jm1]);
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,j]
                     +1.0/dxdy*N_ab[i,jm1]);

else
            # Alternative discretization 

            # -- ux terms -- 
            
            # ux(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(ip1,j,nx,ny);
            Av[k] = ( 4.0/dxdx * N[ip1,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-4.0/dxdx * (N[i,j]+N[ip1,j])
                     -1.0/dydy * (N[i,j]+N[ip1,j])
                     -β_acx[i,j]);

            # ux(i-1,j)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 4.0/dxdx * N[i,j]);

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 1.0/2.0/dydy * (N[i,j]+N[ip1,j])
                     +1.0/8.0/dydy * (N[ip1,jp1]-N[ip1,jm1]+N[i,jp1]-N[i,jm1]) );

            # ux(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jm1,nx,ny);
            Av[k] = ( 1.0/2.0/dydy * (N[i,j]+N[ip1,j])
                     -1.0/8.0/dydy * (N[ip1,jp1]-N[ip1,jm1]+N[i,jp1]-N[i,jm1]) );
            
            # -- uy terms -- 

            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 3.0/2.0/dxdy * (N[i,j]+N[ip1,j])
                     +1.0/dxdy * (N[ip1,j]-N[i,j])
                     +1.0/8.0/dxdy * (N[ip1,jp1]-N[ip1,jm1]+N[i,jp1]-N[i,jm1]) );
            
            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-3.0/2.0/dxdy * (N[i,j]+N[ip1,j])
                     +1.0/dxdy * (N[ip1,j]-N[i,j])
                     -1.0/8.0/dxdy * (N[ip1,jp1]-N[ip1,jm1]+N[i,jp1]-N[i,jm1]) );
            
            # uy(i+1,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,jm1,nx,ny);
            Av[k] = (-3.0/2.0/dxdy * (N[i,j]+N[ip1,j])
                     -1.0/dxdy * (N[ip1,j]-N[i,j])
                     +1.0/8.0/dxdy * (N[ip1,jp1]-N[ip1,jm1]+N[i,jp1]-N[i,jm1]) );
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 3.0/2.0/dxdy * (N[i,j]+N[ip1,j])
                     -1.0/dxdy * (N[ip1,j]-N[i,j])
                     -1.0/8.0/dxdy * (N[ip1,jp1]-N[ip1,jm1]+N[i,jp1]-N[i,jm1]) );
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

if use_sico_discretization
            
            # Standard (sicopolis) discretization
    
            # -- uy terms -- 
            
            # uy(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jp1,nx,ny);
            Av[k] = ( 4.0/dydy*N[i,jp1]);  
            
            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-4.0/dydy * (N[i,jp1]+N[i,j])
                     -1.0/dxdx * (N_ab[i,j]+N_ab[im1,j])
                     -β_acy[i,j]);

            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 4.0/dydy*N[i,j]);    
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 1.0/dxdx*N_ab[i,j]);     
            
            # uy(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(im1,j,nx,ny);
            Av[k] = ( 1.0/dxdx*N_ab[im1,j]);   
            
            # -- ux terms -- 

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,jp1]
                     +1.0/dxdy*N_ab[i,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,j]
                     -1.0/dxdy*N_ab[i,j]); 

            # ux(i-1,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,jp1,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,jp1]
                     -1.0/dxdy*N_ab[im1,j]);
            
            # ux(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,j]
                     +1.0/dxdy*N_ab[im1,j]);
else
            # Alternative discretization 

            # -- uy terms -- 
            
            # uy(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jp1,nx,ny);
            Av[k] = ( 4.0/dydy * N[i,jp1]);

            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-4.0/dydy * (N[i,j]+N[i,jp1])
                     -1.0/dxdx * (N[i,j]+N[i,jp1])
                     -β_acy[i,j]);

            # uy(i,j-1)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 4.0/dydy * N[i,j]);

            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 1.0/2.0/dxdx * (N[i,j]+N[i,jp1])
                     +1.0/8.0/dxdx * (N[ip1,jp1]-N[im1,jp1]+N[ip1,j]-N[im1,j]) );

            # uy(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(im1,j,nx,ny);
            Av[k] = ( 1.0/2.0/dxdx * (N[i,j]+N[i,jp1])
                     -1.0/8.0/dxdx * (N[ip1,jp1]-N[im1,jp1]+N[ip1,j]-N[im1,j]) );
            
            # -- ux terms -- 

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 3.0/2.0/dxdy * (N[i,j]+N[i,jp1])
                     +1.0/dxdy * (N[i,jp1]-N[i,j])
                     +1.0/8.0/dxdy * (N[ip1,jp1]-N[im1,jp1]+N[ip1,j]-N[im1,j]) );
            
            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-3.0/2.0/dxdy * (N[i,j]+N[i,jp1])
                     +1.0/dxdy * (N[i,jp1]-N[i,j])
                     -1.0/8.0/dxdy * (N[ip1,jp1]-N[im1,jp1]+N[ip1,j]-N[im1,j]) );
            
            # ux(i-1,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,jp1,nx,ny);
            Av[k] = (-3.0/2.0/dxdy * (N[i,j]+N[i,jp1])
                     -1.0/dxdy * (N[i,jp1]-N[i,j])
                     +1.0/8.0/dxdy * (N[ip1,jp1]-N[im1,jp1]+N[ip1,j]-N[im1,j]) );
            
            # ux(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 3.0/2.0/dxdy * (N[i,j]+N[i,jp1])
                     -1.0/dxdy * (N[i,jp1]-N[i,j])
                     -1.0/8.0/dxdy * (N[ip1,jp1]-N[im1,jp1]+N[ip1,j]-N[im1,j]) );

end

            # [u] value
            u[nr] = uy[i,j];

            # [b] value 
            b[nr] = taud_acy[i,j];
            
        end
    end

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

    # Define output velocity arrays with new solution
    ux1 = fill(0.0,nx,ny);
    uy1 = fill(0.0,nx,ny);

    for i = 1:nx
        for j = 1:ny
            n = ij2n_ux(i,j,nx,ny);
            ux1[i,j] = unew[n];
            n = ij2n_uy(i,j,nx,ny);
            uy1[i,j] = unew[n];
        end
    end

    return ux1, uy1
end

function calc_picard_convergence_metric(ux,uy,ux0,uy0;mask_acx=nothing,mask_acy=nothing,vel_tol=1e-1,du_reg=1e-5)

    if isnothing(mask_acx)
        # Assume all points are valid
        kkx = findall(abs.(ux) .> vel_tol);
    else
        kkx = findall(abs.(ux) .> vel_tol .&& mask_acx);
    end
    
    if isnothing(mask_acy)
        # Assume all points are valid
        kky = findall(abs.(uy) .> vel_tol);
    else
        kky = findall(abs.(uy) .> vel_tol .&& mask_acy);
    end
    
    res1  = sqrt( sum((ux[kkx] .- ux0[kkx]).^2) + sum((uy[kky] .- uy0[kky]).^2) );
    res2  = sqrt( sum((ux[kkx] .+ ux0[kkx]).^2) + sum((uy[kky] .+ uy0[kky]).^2) );
    resid = 2.0*res1/(res2+du_reg)

    return resid
end

println("Loaded pagos-base.jl.")


