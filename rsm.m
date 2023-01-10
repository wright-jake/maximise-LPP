function [opt_xval, opt_zval] = rsm(A,b,z_coeffs,basic_var,artificial_var)
    %INITIALISE VARIABLES
    %Find A hat and B hat
    A_hat=[A;-z_coeffs'];
    b_hat=[b;0];
    %Take the basic variable columns from A to compute B inverse
    B=A(:,basic_var);
    invB=inv(B);
    %Find the initial BFS and corresponding coefficients
    initial_BFS=zeros(length(z_coeffs),1);
    initial_BFS(basic_var)=invB*b;
    c_B=z_coeffs(basic_var);
    %Compute c_B*invB
    cBBi=(c_B'*invB)';
    %Compute the inverse B hat matrix
    zero=zeros(height(invB),1);
    invB_hat=[invB zero;cBBi' 1];
    %Set degenerate equal to zero as default
    degenerate=0;
    
    %BEGIN ITERATIONS
    %Add a stopping condition for the number of iterations
    max_iterations = 20;        
    for n = 1:max_iterations 
        %Compute z_j-c_j (objective row)
        btm = invB_hat(height(invB_hat),:);
        obj_row=btm*A_hat;
        %Check if current solution is optimal or not, if all elements are
        %positive then it is optimal
        j = find(obj_row < 0);
        if (isempty(j))
            %Return final basic variables
            basic_var
            n;
            %Return the worth
            worth=invB_hat(height(invB_hat),1:(width(invB_hat)-1))
            %Optimal x values
            xB_hat = invB_hat*b_hat;
            opt_xval = xB_hat(1:(height(xB_hat)-1),:);
            %Optimal z value
            opt_zval = xB_hat(height(xB_hat),:);
            %Check for alternate optima
            alt_opt = obj_row;
            alt_opt(basic_var)=[];
            alt_opt=find(alt_opt== 0);
            if (~isempty(alt_opt))
                disp('The LPP has an alternate optima')
                return;
            end
            %Check for degenerate solution
            if degenerate == 1
                disp('The solution to the LPP is degenerate')
                return;
            end
            %Check for infeasible solution
             infeasible = intersect(basic_var, artificial_var);
             if (~isempty(infeasible))
                 opt_xval = NaN;
                 opt_zval = NaN;
                 disp('The solution to the LPP is infeasible')
                 return;
             end
            return;
        end
        %Find the index of the entering variable
        [~,I] = min(obj_row);
        j_in = I(1);
        %Compute x_(j_in)
        xj_hat = invB_hat*A_hat(:,j_in);
        %Check for unbounded solution
        u=xj_hat(xj_hat>0);
        if (isempty(u))
            opt_xval = NaN;
            opt_zval = NaN;
            disp('The LPP has an unbounded solution')
            return;
        end
        %Compute x_B hat
        xB_hat = invB_hat*b_hat;
        %Find minimum from the ratio division
        I = find(xB_hat>0);
        x = xB_hat(I)./xj_hat(I);
        %Find the index of the exiting variable
        j_out=find(x==(min(x(x>0))));
        %Test for degeneracy
        [r,~]=size(j_out);
        if r>1
            degenerate=1
        end
        j_out=j_out(1);
        %Create the pivot row
        invB_hat(j_out,:)=invB_hat(j_out,:)./xj_hat(j_out);
        %Iterate through each row of inverse B hat except the pivot row
        [rows, ~] = size(invB_hat);
        for n=1:rows
            if (n~=j_out)
                %Compute new current row by row operations
                invB_hat(n,:)=(invB_hat(j_out,:)*-(xj_hat(n)))+invB_hat(n,:);
            end
        end
        %Change the basic variables depending on what has left and 
        % what has entered
        basic_var(:,j_out) = j_in;
    end
end

%initial test
% A=[2 -1 2 1 0; 1 0 4 0 1]
% b=[2 4]'
%z_coeffs=[6 -2 3 0 0]'
% basic_var=[4 5]
%artificial_var=[]

%ready mix problem
% A=[6 4 1 0 0 0; 1 2 0 1 0 0; -1 1 0 0 1 0; 0 1 0 0 0 1]
% b=[24 6 1 2]'
% z_coeffs=[5 4 0 0 0 0]'
% basic_var=[3 4 5 6]
% artificial_var=[]

%big M example
%A=[1 2 -1 0 1 0; 4 1 0 -1 0 1]
%b=[7 6]'
%M=100
%z_coeffs=[-1 -1 0 0 -M -M]'
%basic_var=[5 6]
%artificial_var=[5 6]

%coursework test
%A=[5 8 0 2 0 0 0 1 0 0 0; 1 0 1 0 1 1 0 0 1 0 0; 2 7 1 0 0 1 0 0 0 1 0; 1 1 0 0 0 1 -1 0 0 0 1]
%b=[20 41 30 12]'
%M=100
%z_coeffs=[7 2 3 1 1 1 0 0 0 -M -M]'
%basic_var=[8 9 10 11]
%artificial_var=[10 11]

%youtube example q3
%A=[2 5 1 0 0; 1 1 0 -1 1]
%b=[6 2]'
%M=100
%z_coeffs=[1 1 0 0 -M]'
%basic_var=[4 5]
%artificial_var=[5]

%degeneracy example
% A=[1 4 1 0; 1 2 0 1]
% b=[8 4]'
% z_coeffs=[3 9 0 0]'
% basic_var=[3 4]
% artificial_var=[]

%degeneracy example
%A=[30 20 1 0; 5 10 0 1]
%b=[300 110]'
%z_coeffs=[6 8 0 0]'
%basic_var=[3 4]
%artificial_var=[]

%unbounded example
%A=[1 -1 1 0; 2 0 0 1]
%b=[10 40]'
%z_coeffs=[2 1 0 0]'
%basic_var=[3 4]
%artificial_var=[]

%infeasible example
%A=[2 1 1 0 0; 3 4 0 -1 1]
%b=[2 12]'
%M=100
%z_coeffs = [3 2 0 0 -M]'
%basic_var = [3 5]
%artificial_var=[5]

%alternate optima example
%A=[1 2 1 0; 1 1 0 1]
%b=[5 4]'
%z_coeffs = [2 4 0 0]'
%basic_var = [3 4]
%artificial_var=[]

%alternate optima example
%A=[10 5 1 0 0; 6 10 0 1 0; 4 12 0 0 1]
%b=[50 60 48]'
%z_coeffs = [20 10 0 0 0]'
%basic_var = [3 4 5]
%artificial_var=[]

%example question
%A=[1 1 0 1 0 0 0; 1 0 1 0 0 1 0; 0 1 1 0 -1 0 1]
%b=[20 5 10]'
%z_coeffs=[1 -1 3 0 0 -M -M]'
%basic_var = [4 6 7]
%artificial_var = [6 7]

%[opt_xval, opt_zval] = rsm(A,b,z_coeffs,basic_var,artificial_var)

%A=[5 8 0 2 0 0 0 1 0 0 0; 1 0 1 0 1 1 0 0 1 0 0; 2 7 1 0 0 1 0 0 0 1 0; 1 1 0 0 0 1 -1 0 0 0 1]
%b=[20 41 30 12]'
%M=100
%z_coeffs=[7 2 3 1 1 1 0 0 0 -M -M]'
%basic_var=[8 9 10 11]
%artificial_var=[10 11]


