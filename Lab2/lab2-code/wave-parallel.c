/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

/* #define WRITE_TO_FILE */
/* #define VERIFY */
double timer();
double initialize(double x, double y, double t);
void save_solution(double *u, int Ny, int Nx, int n);

int main(int argc, char *argv[])
{
  int N,Nt;
  int nx, ny;
  int startx, starty;
  double dt, dx, lambda_sq;
  double *u;
  double *u_old;
  double *u_new;
  double begin,end;

  int dims[2]={0,0};
  int periods[2]={0,0};
  int reorder=1;
  int ndim=2;
  int nproc;
  int rank;
  int coords[2];
  MPI_Comm proc_grid;

  MPI_Request  req_send[4];
  MPI_Request  req_recv[4];
  MPI_Status   stats[4];

  MPI_Datatype column_type;
  int p_left,p_right,p_up,p_down;

#ifdef WRITE_TO_FILE
  double *uglob=NULL;
  MPI_Request  request;
  MPI_Status   status;
  MPI_Datatype block_type; 	/* custom data type for sending my block */
  int          *nys=NULL;
  int          *nxs=NULL;
  MPI_Datatype *block_types=NULL;
#endif


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Dims_create(nproc,ndim,dims);
  MPI_Cart_create(MPI_COMM_WORLD,ndim,dims,periods,reorder,&proc_grid);
  MPI_Comm_rank(proc_grid,&rank);
  MPI_Cart_coords(proc_grid,rank,ndim,coords);

  /* Global sizes */
  N=240;
  if(argc>1)
    N = atoi(argv[1]);
  dx=1.0/(N-1);
  dt=0.50*dx;
  Nt=1.0/dt;
  lambda_sq = (dt/dx)*(dt/dx);

  if(rank==0) {
    printf("%d x %d matrix\n",N,N);
    printf("%d x %d processes\n",dims[1],dims[0]);
  }


  int q = (N-2)/dims[0];
  int r = (N-2)%dims[0];

  if (r>coords[0]) {
    nx=q+1;
    startx=nx*coords[0];
  }
  else {
    nx=q;
    startx = r*(q+1) + (coords[0]-r)*q;
  }

  q = (N-2)/dims[1];
  r = (N-2)%dims[1];

  if (r>coords[1]) {
    ny=q+1;
    starty=ny*coords[1];
  }
  else {
    ny=q;
    starty = r*(q+1) + (coords[1]-r)*q;
  }

  size_t locsize = (nx+2)*(ny+2)*sizeof(double);
  u     = malloc(locsize);
  u_old = malloc(locsize);
  u_new = malloc(locsize);

  MPI_Type_vector(ny,1,nx+2,MPI_DOUBLE,&column_type);
  MPI_Type_commit(&column_type);

  MPI_Cart_shift(proc_grid,0,1,&p_left,&p_right);
  MPI_Cart_shift(proc_grid,1,1,&p_down,&p_up);

#ifdef WRITE_TO_FILE
	MPI_Type_vector(ny, nx, nx+2, MPI_DOUBLE,
                  &block_type);
	MPI_Type_commit(&block_type);

  if(rank==0) {
    nys     = malloc(nproc*sizeof(int));
    nxs     = malloc(nproc*sizeof(int));
    block_types = malloc(nproc*sizeof(MPI_Datatype));

    for(int i=0; i<dims[1]; i++){
      for(int j=0; j<dims[0]; j++){

        nxs[i*dims[0]+j] = (N-2)/dims[0];
        if (j<(N-2)%dims[0])
          nxs[i*dims[0]+j]++;

        nys[i*dims[0]+j] = (N-2)/dims[1];
        if (i<(N-2)%dims[1])
          nys[i*dims[0]+j]++;

        MPI_Type_vector(nys[i*dims[0]+j], nxs[i*dims[0]+j], N, MPI_DOUBLE,
                        &block_types[i*dims[0]+j]);

        MPI_Type_commit(&block_types[i*dims[0]+j]);
      }
    }
  }

  if(rank==0)
    uglob = malloc(N*N*sizeof(double));
#endif

  /* Setup IC */

  memset(u,0,locsize);
  memset(u_old,0,locsize);
  memset(u_new,0,locsize);

  for(int i = 1; i <= ny; ++i) {
    for(int j = 1; j <= nx; ++j) {
      double x = (j+startx)*dx;
      double y = (i+starty)*dx;
      /* u0 */
      u[i*(nx+2)+j] = initialize(x,y,0);

      /* u1 */
      u_new[i*(nx+2)+j] = initialize(x,y,dt);
    }
  }

#ifdef WRITE_TO_FILE
  if(rank==0)
    memset(uglob,0,N*N*sizeof(double));

  MPI_Isend(&u_new[1*(nx+2)+1],
            1, block_type, 0, 1, proc_grid, &request);

  if(rank==0) {

    int src;
    int src_coords[2];
    int recv_startx, recv_starty;

    recv_starty=0;
    for(int i=0; i<dims[1]; i++){
      recv_startx=0;

      for(int j=0; j<dims[0]; j++){

        src_coords[0]=j;
        src_coords[1]=i;
        MPI_Cart_rank(proc_grid, src_coords, &src);

        MPI_Probe(src, 1, proc_grid, &status);

        MPI_Recv(&uglob[(1+recv_starty)*N+(1+recv_startx)], 1,
                 block_types[i*dims[0]+j], src, 1,
                 proc_grid, MPI_STATUS_IGNORE);

        recv_startx += nxs[i*dims[0]+j];
      }
      recv_starty += nys[i*dims[0]+0];
    }

    save_solution(uglob,N,N,1);
  }

  MPI_Wait(&request,&status);
#endif

#ifdef VERIFY
  double max_error=0.0;
#endif

  /* -------------
       Integrate
     ------------- */

  begin=timer();
  for(int n=2; n<Nt; ++n) {

    /* swap ptrs */
    double *tmp = u_old;
    u_old = u;
    u = u_new;
    u_new = tmp;

  /* communicate u which was calculated previous step */
    int ireq=0;

    /* to/from right */
    if(coords[0] < (dims[0]-1)) {
      MPI_Isend(&u[1*(nx+2)+nx],1,column_type,p_right,n,proc_grid,&req_send[ireq]);
      MPI_Irecv(&u[1*(nx+2)+nx+1],1,column_type,p_right,n,proc_grid,&req_recv[ireq]);
      ireq++;
    }

    /* to/from left */
    if(coords[0] > 0) {
      MPI_Isend(&u[1*(nx+2)+1],1,column_type,p_left,n,proc_grid,&req_send[ireq]);
      MPI_Irecv(&u[1*(nx+2)+0],1,column_type,p_left,n,proc_grid,&req_recv[ireq]);
      ireq++;
    }

    /* to/from up */
    if(coords[1] < (dims[1]-1)) {
      MPI_Isend(&u[ny*(nx+2)+1],nx,MPI_DOUBLE,p_up,n,proc_grid,&req_send[ireq]);
      MPI_Irecv(&u[(ny+1)*(nx+2)+1],nx,MPI_DOUBLE,p_up,n,proc_grid,&req_recv[ireq]);
      ireq++;
    }

    /* to/from down */
    if(coords[1] > 0) {
      MPI_Isend(&u[1*(nx+2)+1],nx,MPI_DOUBLE,p_down,n,proc_grid,&req_send[ireq]);
      MPI_Irecv(&u[1],nx,MPI_DOUBLE,p_down,n,proc_grid,&req_recv[ireq]);
      ireq++;
    }

    MPI_Waitall(ireq,req_recv,stats);
    MPI_Waitall(ireq,req_send,stats);


    for(int i = 1; i <= ny; ++i) {
      for(int j = 1; j <= nx; ++j) {

        u_new[i*(nx+2)+j] = 2*u[i*(nx+2)+j]-u_old[i*(nx+2)+j]
          +lambda_sq* (u[(i+1)*(nx+2)+j] + u[(i-1)*(nx+2)+j]
                       + u[i*(nx+2)+j+1] + u[i*(nx+2)+j-1]
                       -4*u[i*(nx+2)+j]);
      }
    }

#ifdef VERIFY
    double error=0.0;
    for(int i = 1; i <= ny; ++i) {
      for(int j = 1; j <= nx; ++j) {

        double e = fabs(u_new[i*(nx+2)+j]-initialize((j+startx)*dx,(i+starty)*dx,n*dt));
        if(e>error)
          error = e;
      }
    }
    if(error > max_error)
      max_error=error;
#endif

#ifdef WRITE_TO_FILE
    MPI_Isend(&u_new[1*(nx+2)+1],
              1, block_type, 0, n, proc_grid, &request);

    if(rank==0) {

      int src;
      int src_coords[2];
      int recv_startx, recv_starty;

      recv_starty=0;
      for(int i=0; i<dims[1]; i++){
        recv_startx=0;

        for(int j=0; j<dims[0]; j++){

          src_coords[0]=j;
          src_coords[1]=i;
          MPI_Cart_rank(proc_grid, src_coords, &src);

          MPI_Probe(src, n, proc_grid, &status);

          MPI_Recv(&uglob[(1+recv_starty)*N+(1+recv_startx)], 1,
                   block_types[i*dims[0]+j], src, n,
                   proc_grid, MPI_STATUS_IGNORE);

          recv_startx += nxs[i*dims[0]+j];
        }
        recv_starty += nys[i*dims[0]+0];
      }

      save_solution(uglob,N,N,n);
    }

    MPI_Wait(&request,&status);
#endif
  }

  end=timer();

#ifdef VERIFY
  double global_max_error;
  MPI_Reduce(&max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#endif

  if(rank==0) {
    printf("Time elapsed: %g s\n",end-begin);
#ifdef VERIFY
    printf("Maximum error: %g\n",global_max_error);
#endif
  }

#ifdef WRITE_TO_FILE
  if(rank==0) {
    free(uglob);
    free(nys);
    free(nxs);
    free(block_types);
  }
#endif

  free(u);
  free(u_old);
  free(u_new);

  MPI_Finalize();

  return 0;
}

double timer()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

double initialize(double x, double y, double t)
{
  double value = 0;
#ifdef VERIFY
  /* standing wave */
  value=sin(3*M_PI*x)*sin(4*M_PI*y)*cos(5*M_PI*t);
#else
  /* squared-cosine hump */
  const double width=0.1;

  double centerx = 0.25;
  double centery = 0.5;

  double dist = sqrt((x-centerx)*(x-centerx) +
                     (y-centery)*(y-centery));
  if(dist < width) {
    double cs = cos(M_PI*dist/width);
    value = cs*cs;
  }
#endif
  return value;
}

void save_solution(double *u, int Ny, int Nx, int n)
{
  char fname[50];
  sprintf(fname,"solution-%d.dat",n);
  FILE *fp = fopen(fname,"w");

  fprintf(fp,"%d %d\n",Nx,Ny);

  for(int j = 0; j < Ny; ++j) {
    for(int k = 0; k < Nx; ++k) {
      fprintf(fp,"%e\n",u[j*Nx+k]);
    }
  }

  fclose(fp);
}
