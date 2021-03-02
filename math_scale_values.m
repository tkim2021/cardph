function [ newValue ] = math_scale_values( originalValue, minOriginalRange, maxOriginalRange, minNewRange, maxNewRange )
%   MATH_SCALE_VALUES
%   Converts a value from one range into another
%       (maxNewRange - minNewRange)(originalValue - minOriginalRange)
%    y = ----------------------------------------------------------- + minNewRange      
%               (maxOriginalRange - minOriginalRange)
newValue = minNewRange + (((maxNewRange - minNewRange) * (originalValue - minOriginalRange))/(maxOriginalRange - minOriginalRange));
end