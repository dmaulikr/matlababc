function [Et,Hr,Hz] = OUTbc(Et,Hr,Hz,ratio)

Et(:,end) = ratio*Et(:,end-1);
Hr(:,end) = ratio*Hr(:,end-1);
Hz(:,end-1:end) = 0.0;