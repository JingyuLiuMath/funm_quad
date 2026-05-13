function [w,H,h] = arnoldi_last_hm_orth_update(m, w, H, h, AV_big)

global V_big;

% Note:
% AV_big = [V_big(:,1:m),w]*[H;h*em']
% AV_big' * V_big(:, 1 : m) = [H;h*em']'*[V_big(:,1:m),w]'*V_big(:,1:m)
%
% Let 
%      [Q,R] = qr([V_big(:,1:m),w],0);
% then
% AV_big' * V_big(:, 1 : m) = [H;h*em']'*[V_big(:,1:m),w]'*V_big(:,1:m)
%                           = [H;h*em']'*R'*Q'*Q(:,1:m)*R(1:m,1:m)
%                           = [H;h*em']'*R'*eye(m+1,m)*R(1:m,1:m)
%                           = [H;h*em']'*R'*R(1:m+1,1:m)
%
% AV_big'*w = [H;h*em']'*[V_big(:,1:m),w]'*w
%           = [H;h*em']'*R'*Q'*w
%           = [H;h*em']'*R'*R(:,end)

em = zeros(m,1); em(m)=1;
[~,R] = qr([V_big(:,1:m),w],0);

%c = ([H;h*em']'*R'*R(1:m+1,1:m)) \ ([H;h*em']'*R'*R(:,end));

% We can also avoid forming R'*R by partitioning R and simplifying:
R11 = R(1:m,1:m);
r   = R(1:m,end);
rho = R(end,end);
M = H'*R11' + h*(em*r');
c = R11\(r +  M \ (h*rho^2*em));

% full formula (too expensive)
%c = (AV_big' * V_big(:, 1 : m)) \ (AV_big' * w); 

w = w - V_big(:, 1 : m) * c;
H(1:m, m) = H(1:m, m) + c * h;

norm_w = norm(w);
w = w / norm_w;
h = h * norm_w;

end