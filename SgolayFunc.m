
function Input_dynamics = SgolayFunc(Input_dynamics0)

order = 20;
framelen = 51; % should be an odd number

b_fi = sgolay(order,framelen);

for ch = 1:114
    InputCh = Input_dynamics0(:,ch);
    ycenter = conv(InputCh,b_fi((framelen+1)/2,:),'valid');
    
    ybegin = b_fi(end:-1:(framelen+3)/2,:) * InputCh(framelen:-1:1);
    yend = b_fi((framelen-1)/2:-1:1,:) * InputCh(end:-1:end-(framelen-1));
    
    Input_dynamics(:,ch) = [ybegin; ycenter; yend];
    
end

end
