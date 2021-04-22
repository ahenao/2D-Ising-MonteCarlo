#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

/*param filename :the name of the file to open
param rows       :number of rows in the file data
param cols       :number of columns in the files data
return float**   :a 2d float array holding the read in values
			or NULL on failure.*/
float** read_file(std::string filename,int rows,int cols)
{
	std::fstream file;//create a stream for the file
	file.open(filename.c_str(), std::ios::in);//open the file to read in
	if(!file.is_open()){return 0;}//if the file failed to open return NULL

        //float** is just a 2d array like float values[cols][rows]
	//create the column pointers
	float** floats = new float*[cols+1];
	
	//create the row pointers
	for(int i = 0; i <cols;++i){ floats[i] = new float[rows+1]; }

	//read each through each row
	for(int i =0;i<rows;++i)
	{
		//read the values in this row and push into the correct
		//column.floats is [cols][rows]
		for(int j =0;j<cols;++j)//push into the col
		{ file >>floats[j][i]; }
	}
	file.close();//close the file

	return floats;//return the 2d array
}

float mean_j(float *floats[], int dif)
{
  double prom=0,count=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      if(i!=dif){
	prom+=floats[2][i];
	count++;}
       
    }
  mu=prom/count++;
  return mu;
}


float meanabs_j(float *floats[], int dif)
{
  double prom=0,count=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      if(i!=dif){
	if(floats[2][i]>0){
	  prom+=floats[2][i];}
	else{
	  prom-=floats[2][i];}
	count++;}
       
    }
  mu=prom/count++;
  return mu;
}

float mean2_j(float *floats[], int dif)
{
  double prom=0,count=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      if(i!=dif){
	prom+=floats[2][i]*floats[2][i];
	count++;}
       
    }
  mu=prom/count++;
  return mu;
}

float mean3_j(float *floats[], int dif)
{
  double prom=0,count=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      if(i!=dif){
	prom+=floats[2][i]*floats[2][i]*floats[2][i];
	count++;}
       
    }
  mu=prom/count++;
  return mu;
}


float mean4_j(float *floats[], int dif)
{
  double prom=0,count=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      if(i!=dif){
	prom+=floats[2][i]*floats[2][i]*floats[2][i]*floats[2][i];
	count++;}
       
    }
  mu=prom/count++;
  return mu;
}

float stdv_j(float *floats[], int dif)
{
  double prom=0,sum2=0,count=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      prom+=floats[2][i];
      sum2+=floats[2][i]*floats[2][i];
      count++;
    }
  mu=sum2/count-(prom*prom)/(count*count);
  // mu=mu/(count-1);
  return mu;
}


float binder_j(float *floats[], int dif)
{
  double summ4=0,summ2=0,count=0,a=0;
  float mu=0;
  for(int i = int(1e3); i <= 1e5; i++)
    {
      a=floats[2][i]*floats[2][i];
      summ2+=a;
      summ4+=a*a;
      count++;
    }
  summ4=summ4/count;
  summ2=summ2/count;
  mu=0.5*(3-((summ4)/(summ2*summ2)));
  return mu;
}



float variance_jack(float *thetai, double prom)
{
  double suma2=0, count=0,a=0;
  float mu=0;
    for(int p = 0; p <= 1e5-1e3; p++)
    {
      a=thetai[p]-prom;
      suma2+=a*a;
      count++;
    }
    mu=(count-1)*suma2/count;
    return mu;
}



int main()
{
  int rows = int(1e5);//number of rows in the file
  int cols = 3;//number of columns in the file
  float** floats;//2d array
  
  /*the func read_file returns null on failure, so set floats to the 
    value returned (which should be a 2d array) and check for failure.
    if the function failed exit main.*/ 
  if( !(floats = read_file("Energy_100a",rows,cols) ) ){return 0;}//Change "Energy_100a.txt by the filename"
  int option;
  std::cout<<"Choose the Quantity:\n";
  std::cout<<"\t [1]. < X >\n";
  std::cout<<"\t [2]. < |X| >\n";
  std::cout<<"\t [3]. < X**2 >\n";
  std::cout<<"\t [4]. < X**3 >\n";
  std::cout<<"\t [5]. < X**4 >\n";
  std::cout<<"\t [6]. < X**2 > - < X >**2\n";
  std::cout<<"\t [7]. Binder: 0.5(3 - < X**4 >/< X**2 >**2)\n";
  std::cin>>option;
  
  float mean=0;
  float meani[int(1e5-1e3+1)];
  int k=int(1e3);
  double theta=0;
  double cont=0;
  float real=0;
  float bias=0;
  float var=0;

  if(option==1){
    mean=mean_j(floats,0);}
  else if(option==2){
    mean=meanabs_j(floats,0);}
  else if(option==3){
    mean=mean2_j(floats,0);}
  else if(option==4){
    mean=mean3_j(floats,0);}
  else if(option==5){
    mean=mean4_j(floats,0);}
  else if(option==6){
    mean=stdv_j(floats,0);}
  else if(option==7){
    mean=binder_j(floats,0);}

  
  for(int i = int(1e3); i <=1e5; i++)
    {
      if(option==1){
	  meani[i-k]=mean_j(floats,i);}
	else if(option==2){
	  meani[i-k]=meanabs_j(floats,i);}
	else if(option==3){
	  meani[i-k]=mean2_j(floats,i);}
	else if(option==4){
	  meani[i-k]=mean3_j(floats,i);}
	else if(option==5){
	  meani[i-k]=mean4_j(floats,i);}
	else if(option==6){
	  meani[i-k]=stdv_j(floats,i);}
	else if(option==7){
	  meani[i-k]=binder_j(floats,i);}
      theta+=meani[i-k];
      cont++;
    }
  theta=(theta/cont);

  var=variance_jack(meani, theta);
    
  real=cont*mean-((cont-1)*theta);
  bias=(cont-1)*(theta-mean);
  if(option==1){
    std::cout<<"<X>   =\t"<<real<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}
  if(option==2){
    std::cout<<"<|X|>   =\t"<<real<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}
  if(option==3){
    std::cout<<"<X**2>   =\t"<<real<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}
  if(option==4){
    std::cout<<"<X**3>   =\t"<<real<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}
  if(option==5){
    std::cout<<"<X**4>   =\t"<<real<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}
  if(option==6){
    std::cout<<"<X**2> - < X >**2   =\t"<<real<<"\t  sqrt(stdv)= "<<sqrt(real)<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}
  if(option==7){
    std::cout<<"Binder   =\t"<<real<<"\t"<<var<<"\t"<<sqrt(var)<<"\t"<<bias<<"\n";}


//   float meanabs=0;
//   float meanabsi[int(1e5-1e3+1)];
//   double thetabs=0;
//   cont=0;
//   real=0;
//   bias=0;
//   meanabs=meanabs_j(floats,0);
//   for(int i = int(1e3); i <=1e5; i++)
//     {
//       meanabsi[i-k]=meanabs_j(floats,i);
//       thetabs+=meanabsi[i-k];
//       cont++;
//     }
//   thetabs=(thetabs/cont);

//   real=cont*meanabs-((cont-1)*thetabs);
//   bias=(cont-1)*(thetabs-meanabs);
//   std::cout<<"<|X|> =\t"<<real<<"\t\t"<<bias<<"\n";




//   float mean2=0;
//   float mean2i[int(1e5-1e3+1)];
//   double theta2=0;
//   cont=0;
//   real=0;
//   bias=0;
//   mean2=mean2_j(floats,0);
//   for(int i = int(1e3); i <=1e5; i++)
//     {
//       mean2i[i-k]=mean2_j(floats,i);
//       theta2+=mean2i[i-k];
//       cont++;
//     }
//   theta2=(theta2/cont);
  
//   real=cont*mean2-((cont-1)*theta2);
//   bias=(cont-1)*(theta2-mean2);
//   std::cout<<"<X**2> =\t"<<real<<"\t\t"<<bias<<"\n";

//   float mean3=0;
//   float mean3i[int(1e5-1e3+1)];
//   double theta3=0;
//   cont=0;
//   real=0;
//   bias=0;
//   mean3=mean3_j(floats,0);
//   for(int i = int(1e3); i <=1e5; i++)
//     {
//       mean3i[i-k]=mean3_j(floats,i);
//       theta3+=mean3i[i-k];
//       cont++;
//     }
//   theta3=(theta3/cont);
  
//   real=cont*mean3-((cont-1)*theta3);
//   bias=(cont-1)*(theta3-mean3);
//   std::cout<<"<X**3> =\t"<<real<<"\t\t"<<bias<<"\n";
  
  
//   float mean4=0;
//   float mean4i[int(1e5-1e3+1)];
//   double theta4=0;
//   cont=0;
//   real=0;
//   bias=0;
//   mean4=mean4_j(floats,0);
//   for(int i = int(1e3); i <=1e5; i++)
//     {
//       mean4i[i-k]=mean4_j(floats,i);
//       theta4+=mean4i[i-k];
//       cont++;
//     }
//   theta4=(theta4/cont);
  
//   real=cont*mean4-((cont-1)*theta4);
//   bias=(cont-1)*(theta4-mean4);
//   std::cout<<"<X**4> =\t"<<real<<"\t\t"<<bias<<"\n";


//   float stdv=0;
//   float stdvi[int(1e5-1e3+1)];
//   double thetastdv=0;
//   cont=0;
//   real=0;
//   bias=0;
//   stdv=stdv_j(floats,0);
//   for(int i = int(1e3); i <=1e5; i++)
//     {
//       stdvi[i-k]=stdv_j(floats,i);
//       thetastdv+=stdvi[i-k];
//       cont++;
//     }
//   thetastdv=(thetastdv/cont);
  
//   real=cont*stdv-((cont-1)*thetastdv);
//   bias=(cont-1)*(thetastdv-stdv);
//   std::cout<<"<X**2> - <X>**2 =\t"<<real<<"\t\t"<<bias<<"\n";


//   float binder=0;
//   float binderi[int(1e5-1e3+1)];
//   double theta_binder=0;
//   cont=0;
//   real=0;
//   bias=0;
//   binder=binder_j(floats,0);
//   for(int i = int(1e3); i <=1e5; i++)
//     {
//       binderi[i-k]=binder_j(floats,i);
//       theta_binder+=binderi[i-k];
//       cont++;
//     }
//   theta_binder=(theta_binder/cont);
  
//   real=cont*binder-((cont-1)*theta_binder);
//   bias=(cont-1)*(theta_binder-binder);
//   std::cout<<"Binder =\t"<<real<<"\t\t"<<bias<<"\n";

  return 0;
}
