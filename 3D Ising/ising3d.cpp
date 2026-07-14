#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>
using namespace std;

// Grid Size
const int N = 30;
// Number of Monte Carlo steps
const int steps = 100000;

int lattice[N][N][N];

// Periodic boundary condition
int index(int i) { return (i+N)%N; }

// Creates a random spin lattice
void randomlattice()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0; k<N; k++ )
            {
                 lattice[i][j][k] = (rand()%2)*2 -1;
            }
           
        }
    }
}

// Calculates the magnetisation
double magnetisation()
{
    int sum = 0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {   
             for(int k=0; k<N; k++ )
            {
              sum += lattice[i][j][k];  
            }
            
        }
    }
    return abs((double)sum/(N*N*N));
}

void montecarlo(double T)
{
    for(int t=0;t<N*N*N;t++)
    {
        int i = rand()%N;
        int j = rand()%N;
        int k = rand()%N;

        int s = lattice[i][j][k];
        int nb = lattice[index(i-1)][j][k]+
                 lattice[index(i+1)][j][k]+
                 lattice[i][index(j-1)][k]+
                 lattice[i][index(j+1)][k]+
                 lattice[i][j][index(k+1)]+
                 lattice[i][j][index(k-1)];
                 

        int dE = 2*s*nb;

        if(dE<=0 || (double)rand()/RAND_MAX < exp(-dE/T))   // Accept flip according to Metropolis criterion
        {
            lattice[i][j][k] = -s;
        }
    }
}

void showlattice(int k)
{   
    cout<<"slice z = "<<k<<endl;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(lattice[i][j][k]==1)
            {
                cout << "*";
            }
            else
            {
                cout << ".";
            }
        }
        cout << endl;
    }
}

int bfs(int startx, int starty,int startz,bool visited[N][N][N])   // Breadth First Search used to find connected spin clusters
{
    int qx[N*N*N];
    int qy[N*N*N];
    int qz[N*N*N];
    int front = 0;
    int rear = 0;
    qx[rear] = startx;
    qy[rear] = starty;
    qz[rear] = startz;
    int spintype = lattice[startx][starty][startz];

    rear++;
    visited[startx][starty][startz] = true;
    int clustersize = 0;

    while(front<rear)
    {
        int x = qx[front];
        int y = qy[front];
        int z = qz[front];
        front++;
        clustersize++;

        int nx,ny,nz;

        nx = index(x-1);
        ny = y;
        nz = z;
        if(!visited[nx][ny][nz] && spintype==lattice[nx][ny][nz])
        {
            visited[nx][ny][nz] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            qz[rear] = nz;
            rear++;
        }

        nx = index(x+1);
        ny = y;
        nz = z;
        if(!visited[nx][ny][nz] && spintype==lattice[nx][ny][nz])
        {
            visited[nx][ny][nz] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            qz[rear] = nz;
            rear++;
        }

        nx = x;
        ny = index(y-1);
        nz = z;
        if(!visited[nx][ny][nz] && spintype==lattice[nx][ny][nz])
        {
            visited[nx][ny][nz] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            qz[rear] = nz;
            rear++;
        }

        nx = x;
        ny = index(y+1);
        nz = z;
        if(!visited[nx][ny][nz] && spintype==lattice[nx][ny][nz])
        {
            visited[nx][ny][nz] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            qz[rear] = nz;
            rear++;
        }
        
        nx = x;
        ny = y;
        nz = index(z+1);
        if(!visited[nx][ny][nz] && spintype==lattice[nx][ny][nz])
        {
            visited[nx][ny][nz] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            qz[rear] = nz;
            rear++;
        }
        nx = x;
        ny = y;
        nz = index(z-1);
        if(!visited[nx][ny][nz] && spintype==lattice[nx][ny][nz])
        {
            visited[nx][ny][nz] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            qz[rear] = nz;
            rear++;
        }
    }
    return clustersize;
}

// Calculate cluster size distribution
// Also computes total and largest cluster
void clusteranalysis()
{
    bool visited[N][N][N] = {false};
    int histogram[N*N*N+1] = {0};
    int totalcluster = 0;
    int largestcluster = 0;

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {   
            for(int z=0 ; z<N; z++)
            {

                if(!visited[i][j][z])
                {
                    int size = bfs(i,j,z,visited);
                     histogram[size]++;
                    totalcluster++;
                    if(size>largestcluster)
                    {
                        largestcluster = size;
                    }
                }
            }
        }
    }

    cout<<"Cluster Size Distribution"<<endl;
    for(int s=1;s<=N*N*N;s++)
    {
        if(histogram[s]>0)
        {
            cout<<"size "<<s<<"------>"<<histogram[s]<<" cluster"<<endl;
        }
    }
    cout<<"Total cluster ="<<totalcluster<<endl;
    cout<<"Largest cluster ="<<largestcluster<<endl;
}

// Measures similarity to initial configuration
double timecorrelation(int reference[N][N][N])
{
    double sum = 0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
           for(int k =0 ; k<N; k++)
           {
                sum = sum + reference[i][j][k]*lattice[i][j][k];
           }
         }
    }
    return sum/(N*N*N);
}

// Measures correlation between spins separated by distance r
double spatialcorrelation(int r)
{
    double count = 0;
    double sum = 0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {   
            for(int k =0; k<N; k++)
            {
                sum = sum + lattice[i][index(j+r)][k]*lattice[i][j][k];
                count++;
                sum = sum + lattice[index(i+r)][j][k]*lattice[i][j][k];
                count++;
                 sum = sum + lattice[i][j][index(k+r)]*lattice[i][j][k];
                count++;

            }    
        }
    }
    return sum/count;
}

double structurefactor(int k)
{
    const double PI = 3.141592653589793;
    double Asin = 0;
    double Acos = 0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            
            for(int z= 0 ; z<N; z++)
            {
                double angle = 2*PI*k*(i+j+k)/N;
                Asin = Asin + lattice[i][j][z]*sin(angle);
                Acos = Acos + lattice[i][j][z]*cos(angle);
            }
        }
    }
    return (Acos*Acos + Asin*Asin)/(N*N*N);
}

int main()
{
    srand(time(0));
    int sample =50;
    double temperature[] = {1,4.4,4.45,4.50,4.511,4.55,6.0};
    
    for(int i=0;i<7;i++)  // Simulations below Tc, near Tc and above Tc
    {    
        double M[sample] = {0};
        double T = temperature[i];
        randomlattice();

        cout<<"-------------T="<<T<<"-------------"<<endl;

        for(int j=0;j<steps;j++)
        {
            montecarlo(T);
        }

        cout<<"Magnetisation ="<<magnetisation()<<endl;
        /*showlattice(N/2);*/
        for(int s = 0; s< sample ; s++ )
        {
            for( int z = 0; z< 1000; z++){
                montecarlo(T);

            }
            M[s]=magnetisation();
        }
        double mean =0 ;
        for(int z =0; z< sample;z++)
        {
            mean += M[z];
        }
        mean = mean/ sample;
        double sigma = 0;
        for (int i = 0 ; i<sample; i++)
        {
            sigma += (mean-M[i])*(mean-M[i]);

        }
        sigma= sqrt(sigma/(sample-1));
        double error = sigma/sqrt(sample);
        double sumM2 =0;
        for (int i =0 ; i<sample;i++)
        {
            sumM2 += M[i]*M[i];
        }
        double meanM2 = sumM2/sample;
        double sus = ((N*N*N)/T) *(meanM2-mean*mean);
        cout<<"Susceptibility = "<<sus<<endl;

        cout<<"Average magnitisation = "<<mean<<endl;
        cout<<"Standard deviation = "<<sigma<<endl;
        cout<<"Error Bar = "<<error<<endl;
        cout<<"---------------------------------------"<<endl;

       clusteranalysis();
        
        int reference[N][N][N];
        for(int a=0;a<N;a++)
        {
            for(int b=0;b<N;b++)
            {
                
                for(int c=0; c<N; c++)
                {
                    reference[a][b][c] = lattice[a][b][c];
                }
                
                
            }
        }

        cout<<"---------------------------------------"<<endl;
        cout<<"Time Correlation "<<endl;
        for(int t=0;t<500;t++)
        {
            if(t%50==0)
            {
                cout<<"t ="<<t<<"  ---------> c(t) "<<timecorrelation(reference)<<endl;
            }
            montecarlo(T);
        }
        double c[sample][10];
        for(int k=0; k<sample;k++)
        {

            int reference[N][N][N];
            for(int a=0;a<N;a++)
            {
             for(int b=0;b<N;b++)
                {
                    for(int c=0 ; c<N;c++)
                    {
                           reference[a][b][c] = lattice[a][b][c];
                    }
                 
                 }
            }
            int ti =0;
             for(int t=0;t<500;t++)
            {
                if(t%50==0)
                {
                    c[k][ti]= timecorrelation(reference);
                    ti ++;
                }
                montecarlo(T);
             }

            for(int z=0 ; z <1000;z++)
            {
                montecarlo(T);
            }   
        }
        double meanc[10] ={0};
        double sigmac[10] ={0};
        double errorc[10] ={0};
        for(int k =0 ; k<10; k++)
        {
            for(int s= 0; s< sample ;s++ )
            {
                meanc[k] += c[s][k];
            }
            meanc[k]= meanc[k]/sample;
        }
         for(int t =0 ; t<10; t++)
        {
            for(int s=0; s< sample ;s++ )
            {
                sigmac[t] += (c[s][t]-meanc[t])*(c[s][t]-meanc[t]);
            }
            sigmac[t]= sqrt(sigmac[t]/(sample-1));
        }
        for(int k=0; k<10;k++)
        {
            errorc[k] = sigmac[k]/sqrt(sample);
        }
        for(int z= 0;z<10; z++ )
        {
            cout<<"t = "<<z*50<<"   Mean = "<<meanc[z]<<"    standard deviation = "<<sigmac[z]<<"    Error = "<<errorc[z]<<endl;
        }        
       

        cout<<"---------------------------------------"<<endl;
        cout<<"Spatial Correlation "<<endl;
        for(int r=0;r<10;r++)
        {
            cout<<"r = "<<r<<"------------>  G(r) = "<<spatialcorrelation(r)<<endl;
        }
        double G[sample][10];
        for(int s=0 ; s<sample ; s++)
        {
            for(int i = 0; i<1000;i++)
            {
                montecarlo(T);
            }
            for(int z = 0; z<10;z++)
            {
                G[s][z]= spatialcorrelation(z);
            }
        }
         double meanG[10] ={0};
        double sigmaG[10] ={0};
        double errorG[10] ={0};
        for(int i = 0 ;i<10;i++ )
        {
            for(int k = 0; k<sample;k++)
            {
                meanG[i]=meanG[i]+G[k][i];
            }
            meanG[i]=meanG[i]/sample;
        }
        for(int k =0 ; k<10;k++)
        {
            for(int s= 0 ; s<sample; s++)
            {
                sigmaG[k] += (G[s][k]-meanG[k])* (G[s][k]-meanG[k]);
            }
            sigmaG[k]= sqrt(sigmaG[k]/(sample-1));
        }
        for(int k= 0; k<10; k++)
        {
            errorG[k]= sigmaG[k]/sqrt(sample);
        }
        for(int k=0 ; k<10; k++)
        {
          cout<<"r = "<<k<<"   Mean = "<<meanG[k]<<"    standard deviation = "<<sigmaG[k]<<"    Error = "<<errorG[k]<<endl;  
        }

        cout<<"---------------------------------------"<<endl;
        cout<<"Structure Factor  "<<endl;
        for(int k=0;k<10;k++)
        {
            cout<<"k = "<<k<<"---------------> S(k) = "<<structurefactor(k)<<endl;
        }
        double S[sample][10];
        for( int k =0 ; k<sample; k++)
        {
            for(int z= 0 ; z<1000; z++)
            {
                montecarlo(T);

            }
            for(int i =0 ; i<10; i++)
            {
                S[k][i] = structurefactor(i);
            }
        }
        double meanS[10] ={0};
        double sigmaS[10] ={0};
        double errorS[10] ={0};
        for(int k=0 ; k<10;k++)
        {
            for(int z = 0 ; z<sample; z++)
            {
                meanS[k]= meanS[k]+ S[z][k];
            }
            meanS[k]= meanS[k]/sample;
        }
        for(int i =0 ; i<10 ; i++)
        {
            for(int s =0 ; s<sample; s++)
            {
                sigmaS[i] += (S[s][i]- meanS[i])*(S[s][i]- meanS[i]);
            }
            sigmaS[i]= sqrt(sigmaS[i]/(sample-1));
        }
        for(int j = 0; j<10 ; j++)
        {
            errorS[j] = sigmaS[j]/sqrt(sample);
        }
        for(int k=0 ; k<10; k++)
        {
          cout<<"K = "<<k<<"   Mean = "<<meanS[k]<<"    standard deviation = "<<sigmaS[k]<<"    Error = "<<errorS[k]<<endl;  
        }

    }

    return 0;
}