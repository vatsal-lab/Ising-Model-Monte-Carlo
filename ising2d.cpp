#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>
using namespace std;

// Grid Size
const int N = 50;
// Number of Monte Carlo steps
const int steps = 100000;

int lattice[N][N];

// Periodic boundary condition
int index(int i) { return (i+N)%N; }

// Creates a random spin lattice
void randomlattice()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            lattice[i][j] = (rand()%2)*2 -1;
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
            sum += lattice[i][j];
        }
    }
    return abs((double)sum/(N*N));
}

void montecarlo(double T)
{
    for(int k=0;k<N*N;k++)
    {
        int i = rand()%N;
        int j = rand()%N;

        int s = lattice[i][j];
        int nb = lattice[index(i-1)][j]+
                 lattice[index(i+1)][j]+
                 lattice[i][index(j-1)]+
                 lattice[i][index(j+1)];

        int dE = 2*s*nb;

        if(dE<=0 || (double)rand()/RAND_MAX < exp(-dE/T))   // Accept flip according to Metropolis criterion
        {
            lattice[i][j] = -s;
        }
    }
}

void showlattice()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(lattice[i][j]==1)
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

int bfs(int startx, int starty, bool visited[N][N])   // Breadth First Search used to find connected spin clusters
{
    int qx[N*N];
    int qy[N*N];
    int front = 0;
    int rear = 0;
    qx[rear] = startx;
    qy[rear] = starty;
    int spintype = lattice[startx][starty];

    rear++;
    visited[startx][starty] = true;
    int clustersize = 0;

    while(front<rear)
    {
        int x = qx[front];
        int y = qy[front];
        front++;
        clustersize++;

        int nx,ny;

        nx = index(x-1);
        ny = y;
        if(!visited[nx][ny] && spintype==lattice[nx][ny])
        {
            visited[nx][ny] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            rear++;
        }

        nx = index(x+1);
        ny = y;
        if(!visited[nx][ny] && spintype==lattice[nx][ny])
        {
            visited[nx][ny] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            rear++;
        }

        nx = x;
        ny = index(y-1);
        if(!visited[nx][ny] && spintype==lattice[nx][ny])
        {
            visited[nx][ny] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            rear++;
        }

        nx = x;
        ny = index(y+1);
        if(!visited[nx][ny] && spintype==lattice[nx][ny])
        {
            visited[nx][ny] = true;
            qx[rear] = nx;
            qy[rear] = ny;
            rear++;
        }
    }
    return clustersize;
}

// Calculate cluster size distribution
// Also computes total and largest cluster
void clusteranalysis()
{
    bool visited[N][N] = {false};
    int histogram[N*N+1] = {0};
    int totalcluster = 0;
    int largestcluster = 0;

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(!visited[i][j])
            {
                int size = bfs(i,j,visited);
                histogram[size]++;
                totalcluster++;
                if(size>largestcluster)
                {
                    largestcluster = size;
                }
            }
        }
    }

    cout<<"Cluster Size Distribution"<<endl;
    for(int s=1;s<=N*N;s++)
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
double timecorrelation(int reference[N][N])
{
    double sum = 0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            sum = sum + reference[i][j]*lattice[i][j];
        }
    }
    return sum/(N*N);
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
            sum = sum + lattice[i][index(j+r)]*lattice[i][j];
            count++;
            sum = sum + lattice[index(i+r)][j]*lattice[i][j];
            count++;
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
            double angle = 2*PI*k*(i+j)/N;
            Asin = Asin + lattice[i][j]*sin(angle);
            Acos = Acos + lattice[i][j]*cos(angle);
        }
    }
    return (Acos*Acos + Asin*Asin)/(N*N);
}

int main()
{
    srand(time(0));
    int sample =50;
    double temperature[] = {1.0,1.8,2.2, 2.269,2.3,2.4,3.0,4.0};
    
    for(int i=0;i<8;i++)  // Simulations below Tc, near Tc and above Tc
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
        showlattice();
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
        double M2 =0;
        for(int i =0; i<sample; i++)
        {
            M2+=M[i]*M[i];
        }
        double meanM2 = M2/sample;
        double sus = ((N*N)/T)*(meanM2-mean*mean);

        
        cout<<"Average magnitisation = "<<mean<<endl;
        cout<<"Standard deviation = "<<sigma<<endl;
        cout<<"Error Bar = "<<error<<endl;
        cout<<"Susceptibility = "<<sus<<endl;
        cout<<"---------------------------------------"<<endl;

       clusteranalysis();
        
        int reference[N][N];
        for(int a=0;a<N;a++)
        {
            for(int b=0;b<N;b++)
            {
                reference[a][b] = lattice[a][b];
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

            int reference[N][N];
            for(int a=0;a<N;a++)
            {
             for(int b=0;b<N;b++)
                {
                    reference[a][b] = lattice[a][b];
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