N=500;

xmin=-5;

xmax=5;

uph=1;

downh=0.5;

cfl=0.9;

g=1;

tmax=10;

dt = 0.001;

dx=(xmax)/N;

x=(xmin-dx:dx:xmax+dx);

N =length(x);

for i = 1:N
    
    h(i)= 1+1/(sqrt(2.0*pi))*exp(-(3.5*x(i))^2/2.0);
    
    h_ic(i) = h(i);
    
    u(i) = 0;
end

u_ic = u;

length(N);

plot(x,h,'red');  hold on;

%  u(x <= 2.5) = uph;
%
%  u(x > 5) = downh;

t=0;

hu = h.*u;
    
hu_np1 = hu;

h_np1 = h;

y = hu.^2 + 1/2*g.*h.*h;
    

for Step=1:2000
        
    h_new = h;  hu_new = hu;
    
    hu_f = hu; 
    
    h(1) = h(3);
    
    hu(1) = hu(3);
    
    h(N) = h(N-1);
    
    hu(N) = hu(N-1);
    
    for i=3:N-3
        
        hu_0 = 1/3*hu(i-2) -7/6*hu(i-1) + 11/6*hu(i) ;
        
        hu_1 = -1/6*hu(i-1) +5/6*hu(i) + 1/3*hu(i+1) ;
        
        hu_2 = 1/3*hu(i) +5/6*hu(i+1) - 1/6*hu(i+2) ;
        
        gamma_0 = 1/10; gamma_1 = 3/5; gamma_2 = 3/10;
        
        beta_0 = 13/12*(h(i-2)-2*h(i-1)+h(i))^2 ...
            + 1/4*(h(i-2) - 4*h(i-1) +3*h(i))^2;
        
        beta_1 = 13/12*(h(i-1)-2*h(i)+h(i+1))^2 ...
            + 1/4*(h(i-1) - h(i+1))^2;
        
        beta_2 = 13/12*(h(i)-2*h(i+1) + h(i+2))^2 ...
            + 1/4*(3*h(i) - 4*h(i+1) + h(i+2))^2;
        
        epsilon = 1e-7;
        
        
        w_0_tilde = gamma_0/(epsilon+beta_0)^2;
        w_1_tilde = gamma_1/(epsilon+beta_1)^2;
        w_2_tilde = gamma_2/(epsilon+beta_2)^2;
        
        
        w_0 = w_0_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_1 = w_1_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_2 = w_2_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        
        
        hu_f(i+1) = w_0*hu_0 + w_1*hu_1 + w_2*hu_2;
        
    end
    
    for i=3:N-3
        
        y_0 = 1/3*y(i+3) -7/6*y(i+2) + 11/6*y(i+1) ;
        
        y_1 = -1/6*y(i+2) +5/6*y(i+1) + 1/3*y(i) ;
        
        y_2 = 1/3*y(i+1) +5/6*y(i) - 1/6*y(i-1) ;
        
        
        gamma_0 = 1/10; gamma_1 = 3/5; gamma_2 = 3/10;
        
        
        beta_0 = 13/12*(hu_f(i+3)-2*hu_f(i+2)+hu_f(i+1))^2 ...
            + 1/4*(hu_f(i+3) - 4*hu_f(i+2) +3*hu_f(i+1))^2;
        
        beta_1 = 13/12*(hu_f(i+2)-2*hu_f(i+1)+hu_f(i))^2 ...
            + 1/4*(hu_f(i+2) - hu_f(i))^2;
        
        beta_2 = 13/12*(hu_f(i+1)-2*hu_f(i)+hu_f(i-1))^2 ...
            + 1/4*(3*hu_f(i+1) - 4*hu_f(i) +hu_f(i-1))^2;
        
        epsilon = 1e-7;
        
        
        w_0_tilde = gamma_0/(epsilon+beta_0)^2;
        w_1_tilde = gamma_1/(epsilon+beta_1)^2;
        w_2_tilde = gamma_2/(epsilon+beta_2)^2;
        
        
        w_0 = w_0_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_1 = w_1_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_2 = w_2_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        
        
        y_f(i+1) = w_0*y_0 + w_1*y_1 + w_2*y_2;
        
    end

    
    for i=4:N-3
        
        h_np1(i) = h_new(i) - dt*(hu_f(i+1)-hu_f(i))/dx;
        
        hu_np1(i) = hu_new(i) - dt*(y_f(i+1)-y_f(i))/dx;
        
    end
    
    h = h_np1;
    
    hu = hu_np1;
    
    y = 0;
    
    for i=1:N
        
        if (h~=0)
            
            y(i) = hu(i)^2/(h(i)) + 1/2*g*h(i)*h(i);
            
        end
        
    end
    
    for i=3:N-3
        
        hu_0 = 1/3*hu(i-2) -7/6*hu(i-1) + 11/6*hu(i) ;
        
        hu_1 = -1/6*hu(i-1) +5/6*hu(i) + 1/3*hu(i+1) ;
        
        hu_2 = 1/3*hu(i) +5/6*hu(i+1) - 1/6*hu(i+2) ;
        
        gamma_0 = 1/10; gamma_1 = 3/5; gamma_2 = 3/10;
        
        beta_0 = 13/12*(h(i-2)-2*h(i-1)+h(i))^2 ...
            + 1/4*(h(i-2) - 4*h(i-1) +3*h(i))^2;
        
        beta_1 = 13/12*(h(i-1)-2*h(i)+h(i+1))^2 ...
            + 1/4*(h(i-1) - h(i+1))^2;
        
        beta_2 = 13/12*(h(i)-2*h(i+1) + h(i+2))^2 ...
            + 1/4*(3*h(i) - 4*h(i+1) + h(i+2))^2;
        
        epsilon = 1e-7;
        
        
        w_0_tilde = gamma_0/(epsilon+beta_0)^2;
        w_1_tilde = gamma_1/(epsilon+beta_1)^2;
        w_2_tilde = gamma_2/(epsilon+beta_2)^2;
        
        
        w_0 = w_0_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_1 = w_1_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_2 = w_2_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        
        
        hu_f(i+1) = w_0*hu_0 + w_1*hu_1 + w_2*hu_2;
        
    end
    
    for i=3:N-3
        
        y_0 = 1/3*y(i+3) -7/6*y(i+2) + 11/6*y(i+1) ;
        
        y_1 = -1/6*y(i+2) +5/6*y(i+1) + 1/3*y(i) ;
        
        y_2 = 1/3*y(i+1) +5/6*y(i) - 1/6*y(i-1) ;
        
        
        gamma_0 = 1/10; gamma_1 = 3/5; gamma_2 = 3/10;
        
        
        beta_0 = 13/12*(hu_f(i+3)-2*hu_f(i+2)+hu_f(i+1))^2 ...
            + 1/4*(hu_f(i+3) - 4*hu_f(i+2) +3*hu_f(i+1))^2;
        
        beta_1 = 13/12*(hu_f(i+2)-2*hu_f(i+1)+hu_f(i))^2 ...
            + 1/4*(hu_f(i+2) - hu_f(i))^2;
        
        beta_2 = 13/12*(hu_f(i+1)-2*hu_f(i)+hu_f(i-1))^2 ...
            + 1/4*(3*hu_f(i+1) - 4*hu_f(i) +hu_f(i-1))^2;
        
        epsilon = 1e-7;
        
        
        w_0_tilde = gamma_0/(epsilon+beta_0)^2;
        w_1_tilde = gamma_1/(epsilon+beta_1)^2;
        w_2_tilde = gamma_2/(epsilon+beta_2)^2;
        
        
        w_0 = w_0_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_1 = w_1_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_2 = w_2_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        
        
        y_f(i+1) = w_0*y_0 + w_1*y_1 + w_2*y_2;
        
    end

    
    for i=4:N-3
        
        h_np1(i) = 0.75*h_new(i) + 0.25*(h(i) - dt*(hu_f(i+1)-hu_f(i))/dx);
        
        hu_np1(i) = 0.75*hu_new(i) + 0.25*(hu(i) - dt*(y_f(i+1)-y_f(i))/dx);
        
    end
    
    h = h_np1;
    
    hu = hu_np1;
    
    y = 0;
    
    for i=1:N
        
        if (h~=0)
            
            y(i) = hu(i)^2/(h(i)) + 1/2*g*h(i)*h(i);
            
        end
        
    end
    
    for i=3:N-3
        
        hu_0 = 1/3*hu(i-2) -7/6*hu(i-1) + 11/6*hu(i) ;
        
        hu_1 = -1/6*hu(i-1) +5/6*hu(i) + 1/3*hu(i+1) ;
        
        hu_2 = 1/3*hu(i) +5/6*hu(i+1) - 1/6*hu(i+2) ;
        
        gamma_0 = 1/10; gamma_1 = 3/5; gamma_2 = 3/10;
        
        beta_0 = 13/12*(h(i-2)-2*h(i-1)+h(i))^2 ...
            + 1/4*(h(i-2) - 4*h(i-1) +3*h(i))^2;
        
        beta_1 = 13/12*(h(i-1)-2*h(i)+h(i+1))^2 ...
            + 1/4*(h(i-1) - h(i+1))^2;
        
        beta_2 = 13/12*(h(i)-2*h(i+1) + h(i+2))^2 ...
            + 1/4*(3*h(i) - 4*h(i+1) + h(i+2))^2;
        
        epsilon = 1e-7;
        
        
        w_0_tilde = gamma_0/(epsilon+beta_0)^2;
        w_1_tilde = gamma_1/(epsilon+beta_1)^2;
        w_2_tilde = gamma_2/(epsilon+beta_2)^2;
        
        
        w_0 = w_0_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_1 = w_1_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_2 = w_2_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        
        
        hu_f(i+1) = w_0*hu_0 + w_1*hu_1 + w_2*hu_2;
        
    end
    
    for i=3:N-3
        
        y_0 = 1/3*y(i+3) -7/6*y(i+2) + 11/6*y(i+1) ;
        
        y_1 = -1/6*y(i+2) +5/6*y(i+1) + 1/3*y(i) ;
        
        y_2 = 1/3*y(i+1) +5/6*y(i) - 1/6*y(i-1) ;
        
        
        gamma_0 = 1/10; gamma_1 = 3/5; gamma_2 = 3/10;
        
        
        beta_0 = 13/12*(hu_f(i+3)-2*hu_f(i+2)+hu_f(i+1))^2 ...
            + 1/4*(hu_f(i+3) - 4*hu_f(i+2) +3*hu_f(i+1))^2;
        
        beta_1 = 13/12*(hu_f(i+2)-2*hu_f(i+1)+hu_f(i))^2 ...
            + 1/4*(hu_f(i+2) - hu_f(i))^2;
        
        beta_2 = 13/12*(hu_f(i+1)-2*hu_f(i)+hu_f(i-1))^2 ...
            + 1/4*(3*hu_f(i+1) - 4*hu_f(i) +hu_f(i-1))^2;
        
        epsilon = 1e-7;
        
        
        w_0_tilde = gamma_0/(epsilon+beta_0)^2;
        w_1_tilde = gamma_1/(epsilon+beta_1)^2;
        w_2_tilde = gamma_2/(epsilon+beta_2)^2;
        
        
        w_0 = w_0_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_1 = w_1_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        w_2 = w_2_tilde/(w_0_tilde+w_1_tilde+w_2_tilde);
        
        
        y_f(i+1) = w_0*y_0 + w_1*y_1 + w_2*y_2;
        
    end

    
    for i=4:N-3
        
        h_np1(i) = (h_new(i) + 2*(h(i) - dt*(hu_f(i+1)-hu_f(i))/dx))/3;
        
        hu_np1(i) = (hu_new(i) +2*(hu(i) - dt*(y_f(i+1)-y_f(i))/dx))/3;
        
    end
    
    h = h_np1;
    
    hu = hu_np1;
    
    y = 0;
    
    for i=1:N
        
        if (h(i)~=0)
            
            y(i) = hu(i)^2/(h(i)) + 1/2*g*h(i)*h(i);
            
        end
        
    end
    
    figure(1);
    
    subplot(2,1,1)
    plot(x,h_ic,'red');  hold on;    
    plot(x,h,'blue');  hold off;
    axis([-5 5 0.5 1.5]);

    subplot(2,1,2)
    plot(x,h_ic.*(u_ic),'red'); hold on;
    plot(x,hu,'blue');  hold off;
    axis([-5 5 -0.5 0.5]);     
    
    title(t);
    
    pause (0.01);
    
    for i=1:N
        
        u(i) = 1;
        
        if( h(i)~=0)
            
            u(i) = hu(i)/h(i);
            
        end
        
    end
    
    t = t+dt;
    
end

 figure(1);
    
    subplot(2,1,1)
    plot(x,h_ic,'red');  hold on;    
    plot(x,h,'blue');  hold on;
    axis([-5 5 0.5 1.5]);

    subplot(2,1,2)
    plot(x,h_ic.*(u_ic),'red'); hold on;
    plot(x,hu,'blue');  hold on;
    axis([-5 5 -0.5 0.5]);     
    

     
    subplot(2,1,1)

    data_ic = importdata('ic1.txt');
    
    data_ht = importdata('ht02.txt');
    
    data_hut = importdata('hut02.txt');
    
    plot(data_ic(:,1),data_ic(:,2),'yellow'); hold on;
    
    
    plot(data_ht(:,1),data_ht(:,2),'black'); hold on;
    
    subplot(2,1,2)
    
    plot(data_hut(:,1),data_hut(:,2),'black'); hold on;
