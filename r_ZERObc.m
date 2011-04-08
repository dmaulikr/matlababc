function [Et,Hr,Hz] = r_ZERObc(Et,Hr,Hz,cfl)

Et(:,1:2) = 0.0;
Hr(:,1:2) = 0.0;
Hz(:,2) = Hz(:,2) - 2*cfl*Et(:,2);