
double magPerSite(int arr[],int total_Site){
  double sum=0.0;
  for (int i = 0; i < total_Site; i++) {
    sum=sum+arr[i];
  }
  sum=sum/total_Site;

  return sum;
}



int mod2(int a, int b){
  while (a<0) {
    a=a+b;
  }
  return a%b;
}


void corr_2d(int arr[],int Nx, int Ntau,double corr[]){
    double temp=0.0;
    for (int xC = 0; xC < Nx; xC++) {
      for (int tauC = 0; tauC < Ntau; tauC++) {
          temp=0.0;
          for (int xA = 0; xA < Nx; xA++) {
            for (int tauA = 0; tauA < Ntau; tauA++) {
              temp=temp+arr[mod2(xA+xC,Nx)*Ntau+mod2(tauA+tauC,Ntau)]*arr[xA*Ntau+tauA];
            }
          }
          corr[xC*Ntau+tauC]=temp/(Ntau*Nx);

      }
    }
}
