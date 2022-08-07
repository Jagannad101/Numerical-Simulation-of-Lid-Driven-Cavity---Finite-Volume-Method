
clear all
close all
clc
%% Initial Parameters
Nc=32; %No of cells
Nx=Nc+2; %Number of nodes in the X direction - Staggered grid
Ny=Nc+2; %Number of nodes in the Y direction - Staggered grid
n=Nc+2;
N=Nc*Nc;
x = 0 :1/(Nc-1):1;
y = 0:1/(Nc-1):1;
[ X, Y ] = meshgrid (x,y);

dt=5e-4; %time step

dx=1/Nc;
dy=1/Nc;

Re=100; %Reynolds No

U_lid=1; %top lid velocity
L=1;
toler_main = 1e-8; %tolerance for u and v change to reach steady state
toler_p = 1e-5; %tolerance for pressure Poisson's equation
isConverged=0;
itinit =1;
nmax=500000; %Maximum number of iterations
diff = zeros(Nx,Ny); %Diffusion Coefficient


%% Initialize
u=zeros(Nx,Ny) ;
u_tilda=zeros(Nx,Ny) ;
v=zeros(Nx,Ny) ;
v_tilda= zeros(Nx,Ny) ;
p= zeros(Nx,Ny) ;
p_new= zeros(Nx,Ny) ;
u_old=u;
v_old=v;
fp1 = fopen('history_ghost_32_5e-4.dat', 'w');
fprintf(fp1,'TITLE = "Error History "\n' ) ;
fprintf (fp1,'variables=" Iteration ", "time units", "Error_u", "Error_v"\n');
tic

% Start of main loop
for iter=itinit:nmax
    %% Boundary conditions
    %Velocity
    [u,v]=Boundary_Conditions(u,v,U_lid,Nx,Ny);
    %% x_momentum & y_momentum
    [Hx_prev,Hy_prev]=HxHy(u, v, Nx, Ny, Re, dx, dy);
    [Hx, Hy]=HxHy(u_old,v_old,Nx,Ny,G,Re,dx,dy);
    %% Intermediate velocities
    if iter==itinit
        u_tilda(2:Nx-1, 1:Ny-1) = u(2:Nx-1, 1:Ny-1) + dt*Hx(2:Nx-1, 1:Ny-1); %Euler explicit for the first iteration
        v_tilda(1:Nx-1, 2:Ny-1) = v(1:Nx-1 ,2:Ny-1) + dt*Hy(1:Nx-1, 2:Ny-1);
    else
        u_tilda(2:Nx-1,1:Ny-1) = u(2:Nx-1, 1:Ny-1) + dt/2*(3*Hx(2:Nx-1 ,1:Ny-1)-Hx_prev (2:Nx-1, 1:Ny-1)); % Adams Bashforth from the second iteration onwards
        v_tilda (1:Nx-1,2:Ny-1) = v(1:Nx-1, 2:Ny-1) + dt/2*(3*Hy(1:Nx-1 ,2:Ny-1)-Hy_prev(1:Nx-1, 2:Ny-1));
    end
    %% Pressure Poisson equation coefficients
    AS_p = zeros(Nx,Ny);
    AS_p(:,3:Ny-1) = -dx/dy;
    AN_p = zeros(Nx,Ny);
    AN_p(:,2:Ny-2) = -dx/dy;
    AW_p = zeros(Nx,Ny);
    AW_p(3:Nx-1,:) = -dy/dx;
    AE_p = zeros(Nx,Ny) ;
    AE_p(2:Nx-2, :) = -dy/dx ;
    AP_p =(AW_p + AE_p + AS_p + AN_p);
    for i =2:Nx-1
        for j =2:Ny-1
            %Source term using u_tilda
            S(i,j)= ((u_tilda(i+1, j )- u_tilda(i,j)) *dy + (v_tilda(i,j+1) -v_tilda(i,j)) *dx)/dt ;
        end
    end
    %% Gauss Seidel Method for pressure Poisson Equation
    w=1; %Relaxation factor
    res=100;
    while (res>toler_p )
        for i=2:Nx-1
            for j =2:Ny-1
                p_new(i,j) = (1/AP_p(i,j))*(AP_p(i,j)*(1-w)*p(i,j) + w*(AS_p(i,j)*p_new(i,j-1) +AW_p(i,j)*p_new(i-1, j) + AE_p(i,j)*p(i+1, j) + AN_p(i, j) *p(i, j+1) + S(i, j)));
                diff(i,j)=(AE_p(i, j)*p(i+1, j) + AW_p(i, j)*p(i-1, j) + AN_p(i, j)*p(i, j+1) +AS_p(i, j)*p(i, j-1) -AP_p(i, j)*p(i, j) + S(i, j)).^2;
            end
        end
        p=p_new ;
        res=((1/N)*sum(diff( : )))^ 0.5;
    end
    % End of Gauss Seidel
    
    %% Correcting the velocities
    for i =2:Nx-1
        for j =2:Ny-1
            u(i, j)=u_tilda(i, j) -dt*(p(i,j)-p(i-1, j))/dx;
            v(i, j)=v_tilda(i, j) -dt*(p(i,j)-p(i, j-1))/dy;
        end
    end
    err_u=(u-u_old).^2;
    error_u=(sum(err_u( : ))/N)^ 0.5; %L2 norm of U velocity
    err_v=(v-v_old).^2;
    error_v=(sum(err_v ( : ))/N)^ 0.5; %L2 norm of V velocity
    u_old=u;
    v_old=v;
    
    %% Normalize pressure
    p=p-p(5,5);
    t=dt*(iter-1);
    if ( ( mod(iter, 1)==0) || (iter==itinit) )
        fprintf (fp1 , '%d, %e, %e, %e\n', iter, t , error_u , error_v);
        fprintf ('%d, %e, %e, %e\n' , iter , t , error_u , error_v);
    end
    cond = max( error_u , error_v);
    if (cond < toler_main)
        fprintf (fp1 ,'%d, %e, %e, %e\n' , iter, t, error_u, error_v);
        isConverged =1;
        break ;
    end
    
    %End of time iteration loop
end

toc

%% Printing the values
fclose (fp1);
toc
if isConverged == 0
    fprintf (' Solution failed to become steady in %d iterations !!!' ,nmax);
end

if isConverged == 1
    fprintf ('Solution became steady in %d iterations !!!' , iter);
end

%% velocity in cell centers
uc = 0.5*( u( 2 : end-1, 2: end-1) + u( 3 : end, 2: end-1));
vc = 0.5*( v( 2 : end-1, 2: end-1) + v( 2 : end-1, 3: end));

%% Plotting stream-slice

streamslice(X,Y, uc', vc', 2)
axis ( [ 0 ,L, 0 ,L ] )
xlabel( 'X/L ' )
ylabel( 'Y/L ' )
%% Saving the values of x-coordinate, y-coordinate, u-velocity and v-velocity
fp2 = fopen('cavity_ghost_32_5e-4.dat' , 'w' ) ;
fprintf (fp2, 'TITLE = "Cavity Field Data"\n' ) ;
fprintf (fp2, ' variables="x (m)", "y (m)", "u(m/s)", "v(m/s)"\n');
fprintf (fp2, ' zone T="n=%d"\n' , Nc) ;
fprintf (fp2, ' I= %d J= %d\n' ,Nc, Nc) ;
for i =1: length( x )
    for j =1: length( y )
        fprintf (fp2, '%e, %e, %e, %e\n' , x(i) , y(j), uc(i,j) , vc(i,j) ) ;
    end
end
fclose(fp2);

% Function for setting up the Boundary conditions
function[ u, v ] = Boundary_Conditions(u, v, U_lid, Nx, Ny)
%%Bottom wall
%No slip boundary conditions

for i=1:Nx
    u(i,2)=0;
end

for i =1:Nx
    v(i,2)=0;
end

%%Side walls

for j =1:Ny
    u(2,j)=0;
    u(Nx,j)=0;
end

for j =1:Ny
    v(2,j)=0; % 
    v(Nx-1,j)=0;
end

%%Top wall
for i =1:Nx
    u(i,Ny)=U_lid ;
end

for i =1:Nx
    v(i,Ny)=0;
end
end

% Function for the discretized momentum equations
function[Hx,Hy] = HxHy( u, v, Nx, Ny, Re, dx, dy )

for i =2:Nx-1
    for j =2:Ny-1
        rho_u_e(i,j) = (u(i +1,j) + u(i,j) )/2 ;
        rho_u_w(i,j) = (u(i,j) + u(i-1,j) )/2 ;
        rho_v_n(i,j) = (v(i,j+1) + v(i-1,j+1) )/2 ;
        rho_v_s(i,j) = (v(i-1,j) + v (i,j) )/2 ;
    end
end

% x-momentum equation coefficients
for i =2:Nx-1
    for j =2:Ny-1
        AE_u( i , j ) = (rho_u_e(i ,j)/2- 1/(Re*dx ) )*dy ;
        AW_u( i , j ) = (-rho_u_w( i , j ) /2 -1/(Re*dx ) )*dy ;
        AN_u( i , j ) = ( rho_v_n( i , j ) /2 -1/(Re*dy ) )*dx ;
        AS_u( i , j ) = (-rho_v_s( i , j ) /2 - 1/(Re*dy ) )*dx ;
        AP_u( i , j ) = ( rho_u_e( i , j ) /2 + 1/(Re*dx ) ) *dy + (-rho_u_w ( i , j ) /2 +1/(Re*dx ) ) *dy + ( rho_v_n ( i , j ) /2 + 1/(Re*dy ) ) *dx + (-rho_v_s( i , j ) /2 +1/(Re*dy ) ) *dx;
    end
end

    % Correction due to grid spacing at the boundary
    for i =2:Nx-1
        AS_u( i , 2 )=(-rho_v_s( i , j ) /2 - 2/(Re*dy ) ) *dx ;
        AN_u( i ,Ny-1)=( rho_v_n( i , j ) /2 - 2/(Re*dy ) ) *dx ;
        AP_u( i , 2 )=( rho_u_e( i , j ) /2 + 1/(Re*dx ) ) *dy + (-rho_u_w( i , j ) /2 +1/(Re*dx ) ) *dy + ( rho_v_n ( i , j ) /2 + 1/(Re*dy ) ) *dx + (-rho_v_s ( i , j ) /2 +1/(Re*dy ) ) *dx+(1/Re ) *(dy/dx ) ;
        AP_u( i ,Ny-1)=AP_u( i ,Ny-1)+(1/Re )*(dy/dx ) ;
    end
end

for i =3:Nx-1
    for j =2:Ny-1
        Hx( i , j ) = -(AE_u( i , j )*u( i +1, j ) + AW_u( i , j )*u( i-1, j ) + AN_u( i , j )*u( i ,j+1) + AS_u ( i , j ) *u( i , j-1) + AP_u( i , j ) *u( i , j ) ) /( dx*dy ) ;
    end
end

%y momentum equation
for i =2:Nx-1
    for j =2:Ny-1
        AE_v(i,j) = (rho_u_e(i,j)/2 -1/(Re*dx))*dy ;
        AW_v(i,j) = (-rho_u_w(i,j)/2 - 1/(Re*dx))*dy ;
        AN_v(i,j) = (rho_v_n(i,j)/2 -1/(Re*dy))*dx ;
        AS_v(i,j) = (-rho_v_s(i,j)/2 -1/(Re*dy))*dx ;
        AP_v(i,j) = (rho_u_e(i,j)/2 + 1/(Re*dx))*dy + (-rho_u_w (i,j)/2+ 1/(Re*dx)) *dy +(rho_v_n (i,j)/2 + 1/(Re*dy)) *dx + (-rho_v_s (i,j)/2 + 1/(Re*dy))*dx ;
    end
end
    % correction due to grid spacing at the boundary
    for j =2:Ny-1
        AE_v (2 , j )=( rho_u_e ( i , j )/2 -2/(Re*dx ) ) *dy ;
        AW_v (Nx-1, j )=(-rho_u_w ( i , j )/2 -2/(Re*dx ) )*dy ;
        AP_v( 2 , j )=( rho_u_e ( i , j ) /2 + 1/(Re*dx ) ) *dy + (-rho_u_w ( i , j ) /2 +1/(Re*dx ) ) *dy + ( rho_v_n ( i , j ) /2 + 1/(Re*dy ) ) *dx + (-rho_v_s ( i , j ) /2+1/(Re*dy ) ) *dx+(1/Re ) *(dy/dx );
        AP_v (Nx-1, j )=AP_v(Nx-1, j )+(1/Re ) *(dy/dx ) ;
    end
end

for i =2:Nx-1
    for j =3:Ny-1
        Hy( i , j ) = -(AE_v( i , j ) *v ( i +1, j ) + AW_v( i , j ) *v ( i-1, j ) + AN_v( i , j ) *v ( i ,j+1) + AS_v ( i , j ) *v ( i , j-1) + AP_v( i , j ) *v ( i , j ) ) /( dx*dy ) ;
    end
end

end