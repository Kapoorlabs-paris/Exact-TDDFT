
#include<TDSE.h>
#include<wavefunction.h>
#include<fluid.h>
#include<grid.h>
#include<hamop.h>

#define pi 3.1415926
double alphahat;
double frequref;
int main(int argc, char **argv)
{
  int me=0;



  FILE *file_wfdat;
  
  FILE *file_realpot, *file_kohnshamorbital;
  FILE *file_obser, *file_obser_imag;
  FILE *file_reading,*file_readingheliumplus;
 
  
  // *** create some files with appropriate appendices
  
  char string_wfdat[]=       "res/wf_laser_laser_400_2.dat";
   
  char string_realpot[]=        "res/realpot_laser_400_2.dat";
   
  char string_kohnshamorbital[]=       "res/kohnshamorbital_laser_400_2.dat";
  
  char string_obser[]=       "res/observ_laser_400_2.dat";
  char string_obser_imag[]=  "res/observimag_laser_400_2.dat";
  
  char string_reading[]=     "res/wf_singlet_small.dat";
 

   char string_readingheliumplus[]= "res/wf_heliumplus_1500.dat"; 
 
 
 
  
  file_wfdat = fopen(string_wfdat,"w");
  
    
    
   
  file_realpot = fopen(string_realpot,"w");
   
  file_kohnshamorbital = fopen(string_kohnshamorbital,"w");
  file_obser = fopen(string_obser,"w");
  file_obser_imag = fopen(string_obser_imag,"w");
   file_reading = fopen(string_reading,"r");

  

file_readingheliumplus = fopen(string_readingheliumplus,"r");

   long index, rrindex, xindex, yindex,zindex,index2;
 complex<double> imagi(0.0,1.0);
  double deltx=0.2;
double delty=0.2;
double deltz=1;

  long ngpsxbig=1500;
long ngpsybig=1500;

  long ngpsx=1500;
long ngpsy=1500;
    long ngpsz=1;


  // *** declare grid ***
  grid g;
  g.set_dim(16);            // propagation mode is 3 for xy cartesian
  g.set_ngps(ngpsxbig,ngpsybig,ngpsz);    // N_x, N_y, 1 was 100,100,1
  g.set_delt(deltx,delty,deltz);  // delta_x, delta_y, deltz was 0.1,0.1,0.1
  g.set_offs(ngpsxbig/2,ngpsybig/2,0);    // origin  (usually at N_x/2, N_y/2,N_z/2)
  
  // *** declare smaller grid for reading***
  grid g_small;
  g_small.set_dim(16);
  g_small.set_ngps(ngpsx,ngpsy,ngpsz);
  g_small.set_delt(deltx,delty,deltz);
  g_small.set_offs(ngpsx/2,ngpsy/2,0);
  
   




   // *** declare grid for 1D ***
  grid gone;
  gone.set_dim(15);
  gone.set_ngps(ngpsxbig,1,1);
  gone.set_delt(deltx,1,1);
  gone.set_offs(ngpsxbig/2,0,0);
  
  // *** declare grid for 1D ***
  grid gtwo;
  gtwo.set_dim(15);
  gtwo.set_ngps(ngpsxbig,1,1);
  gtwo.set_delt(deltx,1,1);
  gtwo.set_offs(ngpsxbig/2,0,0);  
  
    // *** declare smaler grid for 1D ***
  grid gone_small;
  gone_small.set_dim(15);
  gone_small.set_ngps(ngpsx,1,1);
  gone_small.set_delt(deltx,1,1);
  gone_small.set_offs(ngpsx/2,0,0);
  
   
  
 
  
  
  // *** declare rest of variables ***
  double imag_timestep=0.5;
  double real_timestep=0.2;
  long   no_of_imag_timesteps=0;
  long   no_of_real_timesteps=9300;
  int    obs_output_every=1;
  long   wf_output_every=465; 
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   // don't touch this
  int    box=50;          //was 50
  double masses[]={1.0,1.0};    // don't touch this
  double charge=0.0;    // don't touch this
  double epsilon=0.00001;
  hamop interactionhamil(gone,vecpot_x,vecpot_y,vecpot_z,interactionpotxy,scalarpotytwo,scalarpotztwo,interactionpotxy,imagpotxtwo,imagpotytwo,field,dftpot);
  hamop hamilton(g,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotx,imagpoty,field,dftpot);
 
 
 
 
 hamop hamiltontwo(gone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
 hamop hamiltonone(gone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
  wavefunction wfinidft(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
   wavefunction wfinidftnophase(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());
   wavefunction wfex(g.ngps_x()*g.ngps_y()*g.ngps_z());
   wavefunction wfex2(g.ngps_x()*g.ngps_y()*g.ngps_z());
  
                  wavefunction wfheliumplus(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
                
   wavefunction wfground(g.ngps_x()*g.ngps_y()*g.ngps_z());
   
   
wavefunction wfread(g.ngps_x()*g.ngps_y()*g.ngps_z());
wavefunction wfreadone(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());

wavefunction autowf2(g.ngps_x()*g.ngps_y()*g.ngps_z());


wavefunction wfdftapproxread(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction wfeverythingevendft(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction wfeverythingodddft(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 wavefunction wfplane(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
   wavefunction wfdftapprox(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
    wavefunction wfdftapprox_old(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction lowerimag(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction upper(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction lower(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
    wavefunction lowerimagnophase(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction uppernophase(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction lowernophase(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction potential(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  wavefunction realpot(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
   
   wavefunction wfini(g.ngps_x()*g.ngps_y()*g.ngps_z());
  
  
  wavefunction correlationexdft(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
  
  wavefunction kohnshamdensity(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
  wavefunction kohnshamdensitynophase(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
 
  wavefunction kohnshamorbital(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
  
  wavefunction kohnshamrealdensity(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
   wavefunction kohnshamrealdensitynophase(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
 
 wavefunction correlationdftapprox(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
  wavefunction correlationexdftapprox(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
 wavefunction alpha(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
 wavefunction alphax(gone.ngps_x()*gone.ngps_y()*gone.ngps_z()) ;
double posofinteresty=2;
 double posofinterestx=3;
  //  double posofinterest=2;
  double posofinterestKHx,posofinterestKHy;

  long gpKHx,gplabx,gpKHy,gplaby;

  
long noofalphas=1;
  double deltaalpha=0.02;
  double alphanull=0.1;
  
 long alphacounter;
  
  long nooffrequs=1;
  double deltafrequ=0.005;
  double frequnull=1.35638;
  long frequcounter;


  long nooftargetenergs=2;
  double deltatargetenergs=0.002;
  double targetenergnull=-0.881875;
  long targetenergcounter;
  double targetenergy;

  wavefunction wfdftpot(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 grid g_auto;
  g_auto.set_dim(3);
  g_auto.set_ngps(ngpsx,ngpsy,nooftargetenergs); // <== z runs over the energies of interest 
  g_auto.set_delt(deltx,delty,1.0);
  g_auto.set_offs(ngpsx/2,ngpsy/2,0);
   wavefunction dftpot(g_auto.ngps_x()*g_auto.ngps_y()*g_auto.ngps_z());
  wavefunction autowf(g_auto.ngps_x()*g_auto.ngps_y()*g_auto.ngps_z());
  wavefunction realdftpot(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
  complex<double> timestep;
  double time=0.0;
  
  long counter_i=0;
  long counter_ii=0;
 

  
  
  // initialization
  long outputofinterest=0;
  
  
//wf.init(g,1,2.0,0.0,0.0);
      // ground
    wfread.init(g_small,99,0.1,0.0,0.0,file_reading,outputofinterest);
    wf.nullify();
    wf.regrid(g,g_small,wfread);
    wf*=1.0/sqrt(wf.norm(g));
   wfground=wf;
 
    wfreadone.init(gone_small,99,0.1,0.0,0.0,file_readingheliumplus,outputofinterest);
    wf.nullify();
    wf.regrid(gone,gone_small,wfreadone);
    wf*=1.0/sqrt(wf.norm(gone));

    wfheliumplus=wf;
    
   
 wf=wfground;

   fclose(file_reading);

  cout << "norm wf    : " << wf.norm(g) <<  "\n";

  

  wavefunction staticpot_x(g.ngps_x());
  staticpot_x.calculate_fixed_potential_array_x(g,hamilton,0.0,me);
  
  wavefunction staticpot_y(g.ngps_y());
  staticpot_y.calculate_fixed_potential_array_y(g,hamilton,0.0,me);
  
  
  
  wavefunction staticpot(g.ngps_x()*g.ngps_y()*g.ngps_z());
  staticpot.calculate_fixed_potential_array(g,hamilton,0.0,me);
   

  

  
     wavefunction staticpot_xy(g.ngps_z()*g.ngps_y()*g.ngps_x());
  staticpot_xy.calculate_fixed_potential_array_xy(g,hamilton,0.0,me);
  
  wavefunction staticpotdft(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  staticpotdft.calculate_fixed_potential_array(gtwo,hamiltontwo,0.0,me);

  wavefunction pot(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  pot.nullify();
 wavefunction tddftpot(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 wavefunction tddftpot_old(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 wavefunction tddftpotactual(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());

 complex<double> dftapproxenerg;
 complex<double> complenerg,complenerg0,dftenerg;
  complex<double> groundstatepop, excitedstatepop,excitedstate1dpop,groundstatepopdft;
 complex<double> KHnorm;
  complex<double> complenergex;
  complex<double>correlationint,correlationexint,correlationdftint,correlationexdftint;
  complex<double>correlationdftapproxint,correlationexdftapproxint;
    // ************* imag timeprop
  long ts;
  long no_of_timesteps=no_of_imag_timesteps;
  for (ts=0; ts<no_of_timesteps; ts++)
    {   
                                                                                  
      cout << "Imag: " << ts <<"  "  <<  " wf  "  << real(complenerg0)  <<  endl; 
           
      counter_i++;
      counter_ii++;
      timestep=complex<double>(0.0*real_timestep,-1.0*imag_timestep);  
      time=-imag(timestep*(complex<double>)(ts));
      
      // and now the actual propagation
 
     wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);
       wf*=1.0/sqrt(wf.norm(g));

  
    
      

      
      if (counter_ii==obs_output_every)
      	{
			 complenerg0=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
			 
 	 
	  
	  	  fprintf(file_obser_imag,"%li %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
			  ts,real(complenerg),imag(complenerg),wf.norm(g),wf.non_ionized(g,box),wf.sing_ionized(g,box),wf.doub_ionized(g,box),wf.expect_x(g),wf.expect_y(g));
	  counter_ii=0;
	 
	};
      
    };
  
  
  fclose(file_obser_imag);
  
  
    kohnshamorbital=wf.dens1d(g);
 
  
 
 wfini=wf;
 wfinidft=kohnshamorbital;

  
 
  wf.dump_to_file(g,file_wfdat,dumpingstepwidth);

 

  kohnshamorbital.dump_to_file(gone,file_kohnshamorbital,dumpingstepwidth);


  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
 	 
       
      
//realdftpot=wf.Static_KS(g);  
 
 

    for (alphacounter=0; alphacounter<noofalphas; alphacounter++)
   {
   for (frequcounter=0; frequcounter<nooffrequs; frequcounter++)
   {
    
      counter_i=0;
      counter_ii=0;
      alphahat=alphacounter*deltaalpha+alphanull;
 
      frequref=frequcounter*deltafrequ+frequnull;
   
      wf=wfini;
      kohnshamorbital=wfinidft;
     
  // ************* real timeprop
  
           
  timestep=complex<double>(real_timestep,0.0);
  
  no_of_timesteps=no_of_real_timesteps;
  
  
 



  for (ts=0; ts<no_of_timesteps; ts++)
    {
      counter_i++;
      counter_ii++;
      time=real(timestep*(complex<double>)(ts));
        cout << "Real: " << ts << " alphahat: " << alphahat <<  " frequ: " << frequref << "  "  << " norm wf: " << wf.norm(g) << "  "  <<  " wfdftexact: "  << kohnshamorbital.norm(gone) <<  "\n";
    
if(ts>=no_of_timesteps-1)
    {
      wf=wfini;
      kohnshamorbital=wfinidft;
      
    };


      wf.propagate(timestep,time+0.5*real(timestep),g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);   
 
 
  
  
   
	
	//    Calculating exact DFT stuff
      lowerimag=kohnshamorbital;
      
kohnshamdensity=wf.denskohnsham(g);
	kohnshamrealdensity=kohnshamdensity.denskohnshamcorrector(gone,g,wfheliumplus,wf);
	kohnshamorbital=kohnshamrealdensity.phasekohnshamorbital(gone,g,wf);
 

      upper=kohnshamorbital;
      lower=lowerimag;

      upper.propagatedft(-0.5*timestep,time+0.5*real(timestep),gone,hamiltontwo,me,vecpotflag,0.0*staticpotdft,pot,charge);
      lower.propagatedft(0.5*timestep,time+0.5*real(timestep),gone,hamiltontwo,me,vecpotflag,0.0*staticpotdft,pot,charge);

      
      realpot=(upper/lower).alog(gone,real(timestep));
    
    
	
 

   if (counter_ii==obs_output_every)
      	{
	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);


	
	  	  fprintf(file_obser,"  %.14le  %.14le %.14le %.14le %.14le %.14le  %.14le  %.14le %.14le\n ",
	    (time+real(timestep)),
	    real(complenerg),
		  imag(complenerg),
		  wf.norm(g),

  
 kohnshamorbital.norm(gone),
		
		  wf.expect_x(g),
		  wf.expect_y(g),
		hamilton.vecpot_x(time,me),
			  
		  kohnshamorbital.expect_x(gone)
		 
		  
		  );
	  
	  counter_ii=0;
	};
     

      if (counter_i==wf_output_every)
       	{

	  counter_i=0;
	  wf.dump_to_file(g,file_wfdat,dumpingstepwidth);
 
	  kohnshamorbital.dump_to_file(gone,file_kohnshamorbital,dumpingstepwidth);
	realpot.dump_to_file(gone,file_realpot,dumpingstepwidth);
  
 	};
      	
      
   	  
 

};
 
 
 };
};
  
    
      	

    




    fclose(file_obser);
   
  

    fclose(file_wfdat);
   
    fclose(file_realpot);
    fclose(file_kohnshamorbital);
    cout << me << ": Hasta la vista, ... " << endl;

 

 }



double vecpot_x(double time, int me)
{
  double result=0.0;

  double frequ=frequref;
  double n=400.0;
  double ampl=alphahat*frequ;
  double phi=0.0;
  double dur=n*2.0*M_PI/frequ;
double gap=2.0*M_PI/frequ;
 double ww=0.5*frequ/n;
  

 
   if ((time>0.0)  )
    {
	//result=-0.001;
	result=ampl*sin(ww*time)*sin(ww*time)*sin(frequ*time-phi);
        //result = 0.0;

};



  
      

 




  return result;

}  


double vvvecpot_x(double time, int me)
{

  double result=0.0;

  double frequ= frequref;
  //double frequ2= 0.137;
  double ramping=20.0*2*M_PI/frequ;
  double downramp=20.0*2*M_PI/frequ;
 
  double constperiod=385.0*2*M_PI/frequ;//392.0*2*M_PI/frequ;
  double dur=ramping+constperiod+downramp;
  double ampl = alphahat*frequ;
 
 
	
	

  
 
 
 /*
 
  //double n= dur*frequ/(2*M_PI);

//double gap=2.0*2*M_PI/frequ;

if (time>0)
{
	result=-0.001;

}

*/

//double ramping2=2.0*2*M_PI/frequ2;
 // double downramp2=2.0*2*M_PI/frequ2;
 
  //double constperiod2=104.0*2*M_PI/frequ2;
//double ampl2 = 5.0*alphahat*frequ2;
 
 
  if (time>0)
    {

      if (time<ramping){
	result=-ampl/ramping*time*cos(frequ*time) + ampl/(frequ*ramping)*sin(frequ*time);
	}
      
      if ( (time>=ramping)&&(time<ramping+constperiod))
      {
	result=-ampl*cos(frequ*time);//-ampl*cos(frequ2*time);
      }
     
      if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp)) 
	{
		  result=ampl/downramp*(time-downramp-constperiod-ramping)*cos(frequ*time)-ampl/(frequ*downramp)*sin(frequ*time);
	}
      
      if ( (time >= ramping+constperiod+downramp )) //&& (time<ramping+constperiod+downramp +gap)) 
      {
      result=0.0;
 }
  };
 

    
  return result;
}  



double vvvvvvecpot_x(double time, int me)
{

  double result=0.0;
  return result;
}  


double alpha_y(double time, int me)
{
  double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

   // <-------- put same in vecpot_y !!!
  vecpotampl = alphahat*frequref;


  result=-vecpotampl/frequref*cos(frequref*time);

  return result;

}  


double alpha_x(double time, int me)
{
 double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

  // <-------- put same in vecpot_y !!!
  vecpotampl = alphahat*frequref;


  result=-vecpotampl/frequref*cos(frequref*time);

  return result;

 

}  

double alpha_z(double time, int me)
{
 

  return alpha_y(time,me);

}  


double vecpot_y(double time, int me)
{
 return vecpot_x(time,me);
}  


double vecpot_z(double time, int me)

{

 
 
    

  return 0.0;
 
 // return vecpot_x(time,me);
}  





double scalarpotx(double x, double y, double z, double time, int me)
{
  double eps=1.0;
 double result;
 
 double slope=0.00522;
 double Static_field=0.0141;
 double field_ramp=27;
 double field_const=270;
 result= -2.0/sqrt(x*x+eps);   //x*Static_field-2.0/sqrt(x*x+eps);
 /*
 if (time>0)
    {
		result=x*Static_field-2.0/sqrt(x*x+eps);

if (time<field_ramp){
	result=x*slope*time-2.0/sqrt(x*x+eps);
	}
	if ((time>=field_ramp) && (time<field_ramp+field_const))
	{
		result=x*Static_field-2.0/sqrt(x*x+eps);
	}
	if(time>=field_ramp+field_const) {
		result= -2.0/sqrt(x*x+eps);
	}
	
   //return -2.0/sqrt(x*x+eps);
        //     return 0;
	};
	 double eps2=0.5034;
     return -3.0/sqrt(x*x+eps2*eps2);
    */
	return result;
}

double scalarpoty(double x, double y, double z, double time, int me)
{
  double eps=1.0;
    
double result;
 double slope=0.00522;
 double Static_field=0.00141;
 double field_ramp=27;
 double field_const=270;
 
 result= -2.0/sqrt(y*y+eps);   //y*Static_field-2.0/sqrt(y*y+eps);
 /*
 if (time>0)
    {
		

if (time<field_ramp){
	result=y*slope*time-2.0/sqrt(y*y+eps);
	}
	if ((time>=field_ramp) && (time<field_ramp+field_const))
	{
		result=y*Static_field-2.0/sqrt(y*y+eps);
	}
	if(time>=field_ramp+field_const) {
		result= -2.0/sqrt(y*y+eps);
	}
	
   //return -2.0/sqrt(x*x+eps);
        //     return 0;
	};
 //return -2.0/sqrt(y*y+eps);
      //         return 0;
      */
      return result;
// double eps2=0.5034;
  //   return -3.0/sqrt(y*y+eps2*eps2);
    
}

double scalarpotz(double x, double y, double z, double time, int me)
{
double eps=1.0;

 return 0.0;
 //return -3.0/sqrt(z*z+eps);
}

double interactionpotxy(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  
  return 1.0/sqrt((x-y)*(x-y)+eps);
   
   double eps2=0.5034;
  
  // return 1.0/sqrt((x-y)*(x-y)+eps2*eps2);
   // return 0.0;
}

double interactionpotyz(double x, double y, double z, double time, int me)
{
  double eps=1;
   return 0.0;  
 //return 1.0/sqrt((z-y)*(z-y)+eps);
 
}
double interactionpotxz(double x, double y, double z, double time, int me)
{
  double eps=1;
    return 0.0;  
 //return 1.0/sqrt((x-z)*(x-z)+eps);
 
}
  
double field(double time, int me)
{
  double result=0.0;


  return result;


}  







double imagpotx(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
//    double ampl=0.0; // switch imaginary potential off
   double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
  
      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };
      
}  


double imagpoty(long xindex, long yindex, long zindex, double time, grid g)
{
  double y;
//  double ampl=0.0; // switch imaginary potential off
  
   double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	    *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };
      
} 

double imagpotz(long xindex, long yindex, long zindex, double time, grid g)
{
  double z;
  // double ampl=0.0; // switch imaginary potential off
   double ampl=0.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       z=((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z())
	    *((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z());
       return ampl*z*z*z*z*z*z*z*z;
    }
  else
    {
       return 0.0;
    };
      
} 


// dft is not used
double dftpot(grid g, double x, double y, double z, double time, int me, 
	      const fluid &v_null, const wavefunction &v_eins)
{
  double result;
  result=0.0;
  return result;
}





double scalarpotxtwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;

  return -2.0/sqrt(x*x+eps);
         //
          //  double eps=0.5034;
     //return -3.0/sqrt(x*x+eps2*eps2);
    
           
}

double scalarpotytwo(double x, double y, double z, double time, int me)
{
   double result=0.0;
  return result;

}

double scalarpotztwo(double x, double y, double z, double time, int me)
{
  double result=0.0;
  return result;
}

double interactionpotxytwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  double interactionpotyztwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  double interactionpotxztwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  

double imagpotxtwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
  
    //double ampl=0.0; // switch imaginary potential off
     double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
  
      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };
      
}  


double imagpotytwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double y;
    double ampl=0.0; // switch imaginary potential off
  // double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	    *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };
      
} 

double imagpotztwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double z;
    double ampl=0.0; // switch imaginary potential off
  // double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       z=((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z())
	    *((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z());
       return ampl*z*z*z*z*z*z*z*z;
    }
  else
    {
       return 0.0;
    };
      
} 
double dfthartree(grid gone, double x, double y, double z, double time, int me, 
	      const wavefunction & wf, double eps, long box)
{
  double result=0.0;
  return result;
}


