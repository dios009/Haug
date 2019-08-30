// Differential cross section of electron-positron Bremsstrahlung - Eberhard Huag

double M0 (double k1, double k2, double ro, double tau){
    return
sqrt(pow(ro,2)-4.)*((9.*tau+4)/(4.*k1*k2)-(3.*(pow(ro,2)+3))/(4.*pow(k1,2))-k2*pow(ro,2)/(4.*pow(k1,3))+pow(ro,2)*(pow(k2,2)/pow(k1,3)-2.*k2*tau/pow(k1,2)+k2*tau/(2.*k1)-1./2.)/(k1*(pow(tau,2)-1.))+(k2/k1-1./(tau+1.))/(2.*(tau+1.))-1./(k1+k2)+((k1*k2-2.5*pow(ro,2))/(tau+1.)-pow(ro,2)/((2.*(tau+1))*(tau+1.))-4*pow(ro,2))/((k1+k2)*(k1+k2))+k1*k2*pow(ro,2)*(4.+1./(tau+1.))/pow((k1+k2),4)+(4.*tau/(k1*k2)-4./pow(k1,2)+4*k2/k1+3.)/pow(ro,2)+(pow(ro,2)-4.)*((pow(ro,2)-4.)*(tau/k2-1./k1)-4*k2)/(12.*k1*pow(ro,2)*pow(ro,2)));
}

double MX1 (double k1, double k2, double ro, double tau){
    return
((tau+1.)*(k2-k1)-k2*(k1+k2))*(sqrt(pow(ro,2)-4.)+(k2-tau+1.)*ro*log((ro*(tau-k2-1.)+sqrt(pow(ro,2)-4.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))/(2.*k1))/sqrt((tau-k2)*(tau-k2)+2*k1-1.))/(sqrt((tau-k2)*(tau-k2)+2*k1-1.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))*((2.*tau+3.)/(k2*pow(ro,2))-2./(k1*k2*pow(ro,2))-1/(k1*pow(ro,2))+(3.*tau-1.)/(2.*k1*k2)+k2/(2.*k1*(tau+1.))-1./k1);
}

double MX2 (double k1, double k2, double ro, double tau){
    return
    2.*ro*(ro*sqrt(pow(ro,2)-4.)*(k1*tau-k2)/(k1+k2)+2.*(-k2*tau+k1+pow(tau,2)-1.)*log((pow(tau,2)-1.+(1./2.)*ro*sqrt((pow(ro,2)-4.)*(pow(tau,2)-1.)))/(k1+k2)-tau)/sqrt(pow(tau,2)-1.))/(pow(k1,2)*(pow(tau,2)-1.))*(tau+1.-1./(2.*(tau+1.)));
}

double MX3 (double k1, double k2, double ro, double tau){
    return
    -ro*((tau+1.)*(k2-k1)-k2*(k1+k2))*(2.*log((ro*(tau-k2-1.)+sqrt(pow(ro,2)-4.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))/(2.*k1))/sqrt((tau-k2)*(tau-k2)+2*k1-1.)+(k2-tau+1.)*ro*sqrt(pow(ro,2)-4.)/(2.*pow(k1,2)))/(sqrt((tau-k2)*(tau-k2)+2*k1-1.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))/(2.*k1);
}

double MX4 (double k1, double k2, double ro, double tau){
    return
    -((k2-tau+1.)*sqrt(pow(ro,2)-4.)*(((tau+1.)*(k2-k1)-k2*(k1+k2))*((tau+1.)*(k2-k1)-k2*(k1+k2))-(1./2.)*pow(ro,2)*(2.*k1*k2*tau-pow(k1,2)-pow(k2,2)))+(((k2-tau+1.)*(k2-tau+1.))*((tau+1.)*(k2-k1)-k2*(k1+k2))*((tau+1.)*(k2-k1)-k2*(k1+k2))-2.*pow(k1,2)*(2.*k1*k2*tau-pow(k1,2)-pow(k2,2)))*ro*log((ro*(tau-k2-1.)+sqrt(pow(ro,2)-4.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))/(2.*k1))/sqrt((tau-k2)*(tau-k2)+2*k1-1.))/((sqrt((tau-k2)*(tau-k2)+2*k1-1.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))*sqrt((tau-k2)*(tau-k2)+2*k1-1.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))*(tau+1.-k2)/(2.*k1*k2*pow(ro,2));
}

double MX5 (double k1, double k2, double ro, double tau){
    return
    -(tau*tau)*ro*(ro*sqrt(pow(ro,2)-4.)*(1.-3.*(2.*k1*k2*tau-pow(k1,2)-pow(k2,2))/(pow(k1,2)*(pow(tau,2)-1.)))+((2.*(tau+1.)-k2*pow(ro,2)/k1-3.*(2.*k1*k2*tau-pow(k1,2)-pow(k2,2))*(k1+k2-(1./2.)*(tau-1.)*pow(ro,2))/(pow(k1,2)*(pow(tau,2)-1.)))*2.)*log((pow(tau,2)-1.+(1./2.)*ro*sqrt((pow(ro,2)-4.)*(pow(tau,2)-1.)))/(k1+k2)-tau)/sqrt(pow(tau,2)-1.))/(pow(k1,2)*(pow(tau,2)-1.));
}

double ML1(double k1, double k2, double ro, double tau){
    return
    ro*log((ro*(tau-k2-1.)+sqrt(pow(ro,2)-4.)*sqrt((tau-k2)*(tau-k2)+2*k1-1.))/(2.*k1))*(((3.*pow(k2,2)+2.*k2-4.)/k1+3.*k2+6.)/(2.*(tau+1.))+(5.*tau-11.*k2+1.)/(2.*k1)-(13.*tau+7.)/(2.*k2)+2.*pow(ro,2)/pow(k1,2)-2.*tau*pow(ro,2)/(k1*k2)+(7.*k1+4*k2-10.*tau+2.+(-k2*tau+k2+pow(k2,2)+4.)/k1-(11.*k1*tau+9.*k1+4.*tau-8.)/k2-8./pow(k1,2)+8.*tau/(k1*k2))/(2.*pow(ro,2)))/sqrt((tau-k2)*(tau-k2)+2*k1-1.);
}

double ML2(double k1, double k2, double ro, double tau){
    return
    4.*ro*log(1.+(k2*pow(ro,2)*(pow(ro,2)-4.)+ro*sqrt(pow(ro,2)-4.)*sqrt(k2*(k2*pow(ro,2)*(pow(ro,2)-4.)+8.*k1*(k1+k2))))/(4.*k1*(k1+k2)))*(k1-tau+2.-2.*tau/k2+((tau-5.*k2*(1./4.))*pow(ro,2)+4.*pow(tau,2)-3.*k2*tau+2.*pow(k2,2))/k1+k2*(tau-k1-2.*k2)/(2.*(tau+1.))-(k2*pow(k2,2)-3.*k2+tau)/(k1*(tau+1.))+(k1*k2+3.*k2*tau-k2-4.*tau-4.*pow(tau,2)+2.)/(2.*(k1+k2))-2.*pow(tau,2)*tau/(k1*(k1+k2))+(k1-1.)/((k1+k2)*(tau+1.)))/sqrt(k2*(k2*pow(ro,2)*(pow(ro,2)-4.)+8.*k1*(k1+k2)));
}

double ML3(double k1, double k2, double ro, double tau){
    return
    ro*log((ro+sqrt(pow(ro,2)-4))*(1./2.))*((1.-3.*k2/k1+1/(k1*k2)-8./(k1+k2)+(2.*(-k1*k2+2.*tau))/pow((k1+k2),2)+4.*k1*k2/pow((k1+k2),3)+4.*k1*k2/pow((k1+k2),4))/(tau+1.)+(1.-2./(k1+k2)-2./pow((k1+k2),2))/((tau+1.)*(tau+1.))-6./k1-8./pow(k1,2)-4.*pow(ro,2)/(k1*k2)+5.*pow(ro,2)/pow((k1+k2),2)-((pow(ro,2)+6.)*2.)*k1*k2/pow((k1+k2),4)+2.*(1./k2-k2-2*(tau+1.))/(k1*pow(ro,2)));
}

double ML4(double k1, double k2, double ro, double tau){
    return
    ro*log((pow(tau,2)-1.+(1./2.)*ro*sqrt((pow(ro,2)-4.)*(pow(tau,2)-1.)))/(k1+k2)-tau)*(2.*pow(k1,2)+5.*k1*k2*(1./2.)-tau*(pow(ro,2)-2.)+(1./2.)*k1*pow(ro,2)+2+2*k2*(-2.*k2*tau+pow(k2,2)-pow(ro,2)+2.*pow(tau,2)-1.)/k1+(-3.*k1*k2*tau+k1*k2+4.*tau*pow(tau,2))/(k1+k2)+(pow(k1,2)-2.*k1-1.+2.*k1*k2/(k1+k2)+2.*k2*(k2+1.)/k1)/(tau+1.)+(3.*tau*pow((k1+k2),2)-2.*k1*k2*(tau+1.)+4.*pow(tau,2)-2.*tau-2.-4.*k2*tau/k1)/pow(ro,2)+k2*(k1+k2)*(tau-k2/k1)*(4.*tau*(tau-k2)/pow(k1,2)-1.)/(pow(tau,2)-1.))/(k1*k2*sqrt(pow(tau,2)-1.));
}

  
double tau (double E2){
    return
-E2; 
}


double k2 (double E2, double k, double p2, double x2){
    return
k*(-p2*x2+E2); 
}

double k1 (double k){
    return
-k; 
}

double ro (double E2, double k, double p2, double x2){
    return
sqrt(2.*(k+1-E2-k*(-p2*x2+E2)));
}


/*
double p2 (double E2){
    return
sqrt(E2*E2-1.);
}
*/

double x2min (double E2, double k, double p2){
    return
((E2-1.)*k+E2+1.)/(p2*k); 
} //version corregida del Maple original

/*
double thetamax (double xmin){
    return
acos(xmin);
}
*/

double emin (double k){
    return
(k*k-1.-k*sqrt(k*k-4.*k))/(2.*k+1.);
} 

double emax (double k){
    return
(k*k-1.+k*sqrt(k*k-4.*k))/(2.*k+1.);
}