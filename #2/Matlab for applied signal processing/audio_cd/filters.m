function y = filters(N,D,x,Fs)

% function y = filters(N,D,x,Fs)  simulates analog filtering of the data in
% vector x with the filter N(s)/D(s) described by vectors N and D to create
% the filtered data y. Fs is the sampling frequency of the input (and
% output). Filtering is performed by convolving the input with an estimate
% of the impulse response of the filter, obtained by partial fraction
% expansion. 
% 
% T DUTOIT, 13:49 12/03/2007

[r,p,k]=residue(N,D);
dominant_pole_damping=min(abs(real(p)));

impulse_response_duration=10/dominant_pole_damping;
t=(0:1/Fs:impulse_response_duration);
imp=real(r.'*exp(p*t));

y=conv(x,imp/Fs);


