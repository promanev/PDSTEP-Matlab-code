function result = calcFit(lPos, rPos,pelPos)
% Function to debug fitness function
leftTarg = [0.360931, 2.2843];
rightTarg = [-0.360931, 2.2843];
initPelvisHeight = 3.42645;

result = (-0.1*abs(leftTarg(1) - lPos(1)) + 1) + (-0.1*abs(leftTarg(2) - ...
    lPos(2)) + 1) + (-0.1*abs(rightTarg(1) - rPos(1)) + 1) + ...
    (-0.1*abs(rightTarg(2) - rPos(2)) + 1) + pelPos / initPelvisHeight;
result = result / 5;