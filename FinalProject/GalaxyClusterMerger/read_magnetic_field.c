#include "hdf5.h"
#include "constants.h"
#include "mangle_names.h"
#include <math.h>
#include <stdlib.h>

static hid_t file_id;
static double *xcoord, *ycoord, *zcoord;

void FTOC(open_mag_file)(char filename[MAX_STRING_LENGTH+1], 
			 int *nx,
			 int *ny,
			 int *nz,
			 double *dx,
			 double *dy,
			 double *dz,
			 double *xmin,
			 double *ymin,
			 double *zmin)
{

  herr_t status;
  hid_t dataset, dataspace;
  hsize_t dims[3], maxdims[3];

  int len = 0;
  char* string_index;
  char local_filename[MAX_STRING_LENGTH+1];

  string_index = filename;
  
  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';

  file_id = H5Fopen(local_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = H5Dopen(file_id, "magnetic_vector_potential_z");

  dataspace = H5Dget_space(dataset);

  H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  H5Sclose(dataspace);
  H5Dclose(dataset);

  xcoord = (double *)malloc(dims[0]*sizeof(double));
  ycoord = (double *)malloc(dims[1]*sizeof(double));
  zcoord = (double *)malloc(dims[2]*sizeof(double));

  dataset = H5Dopen(file_id, "/x");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		    H5S_ALL, H5P_DEFAULT, xcoord);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "/y");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                    H5S_ALL, H5P_DEFAULT, ycoord);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "/z");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                    H5S_ALL, H5P_DEFAULT, zcoord);
  H5Dclose(dataset);

  *dx = xcoord[1]-xcoord[0];
  *dy = ycoord[1]-ycoord[0];
  *dz = zcoord[1]-zcoord[0];

  *xmin = xcoord[0]-0.5*(*dx);
  *ymin = ycoord[0]-0.5*(*dy);  
  *zmin = zcoord[0]-0.5*(*dz);
  
  *nx = dims[0];
  *ny = dims[1];
  *nz = dims[2];

  return;

}

void FTOC(read_density_field)(int *ibegin,
			      int *jbegin,
			      int *kbegin,
			      int *iend,
			      int *jend,
			      int *kend,
			      double *Dp)
{

  hid_t dataset, dataspace, memspace, dxfer_template;

  herr_t status;

  hsize_t start[3], stride[3], count[3], dims[3];

  int rank, ierr;

  rank = 3;

  start[0] = *ibegin-1;
  start[1] = *jbegin-1;
  start[2] = *kbegin-1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = *iend-*ibegin+1;
  count[1] = *jend-*jbegin+1;
  count[2] = *kend-*kbegin+1;

  dims[0] = count[0];
  dims[1] = count[1];
  dims[2] = count[2];

  /* Dp */
  dataset = H5Dopen(file_id, "density_perturb");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, Dp);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return;

}


void FTOC(read_magnetic_field)(int *ibegin,
			       int *jbegin,
			       int *kbegin,
			       int *iend,
			       int *jend,
			       int *kend,
			       double *Ax, 
			       double *Ay,
			       double *Az)
{

  hid_t dataset, dataspace, memspace, dxfer_template;

  herr_t status;

  hsize_t start[3], stride[3], count[3], dims[3];

  int rank, ierr;
  
  rank = 3;

  start[0] = *ibegin-1;
  start[1] = *jbegin-1;
  start[2] = *kbegin-1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = *iend-*ibegin+1;
  count[1] = *jend-*jbegin+1;
  count[2] = *kend-*kbegin+1;

  dims[0] = count[0];
  dims[1] = count[1];
  dims[2] = count[2];

  /* Ax */
  dataset = H5Dopen(file_id, "magnetic_vector_potential_x");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
			       stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
  //                             stride, count, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
		   H5P_DEFAULT, Ax);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  /* Ay */
  dataset = H5Dopen(file_id, "magnetic_vector_potential_y");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
			       stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
  //                             stride, count, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
		   H5P_DEFAULT, Ay);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /* Az */
  dataset = H5Dopen(file_id, "magnetic_vector_potential_z");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
			       stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
  //                             stride, count, NULL);  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
		   H5P_DEFAULT, Az);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return;

} 

void FTOC(read_velocity_field)(int *ibegin,
			       int *jbegin,
			       int *kbegin,
			       int *iend,
			       int *jend,
			       int *kend,
			       double *Ax, 
			       double *Ay,
			       double *Az)
{

  hid_t dataset, dataspace, memspace, dxfer_template;

  herr_t status;

  hsize_t start[3], stride[3], count[3], dims[3];

  int rank, ierr;
  
  rank = 3;

  start[0] = *ibegin-1;
  start[1] = *jbegin-1;
  start[2] = *kbegin-1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = *iend-*ibegin+1;
  count[1] = *jend-*jbegin+1;
  count[2] = *kend-*kbegin+1;

  dims[0] = count[0];
  dims[1] = count[1];
  dims[2] = count[2];

  /* Ax */
  dataset = H5Dopen(file_id, "velocity_x");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
			       stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
  //                             stride, count, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
		   H5P_DEFAULT, Ax);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  /* Ay */
  dataset = H5Dopen(file_id, "velocity_y");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
			       stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
  //                             stride, count, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
		   H5P_DEFAULT, Ay);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /* Az */
  dataset = H5Dopen(file_id, "velocity_z");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
			       stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  //status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start,
  //                             stride, count, NULL);  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
		   H5P_DEFAULT, Az);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return;

} 

void FTOC(close_mag_file)()
{

  free(xcoord);
  free(ycoord);
  free(zcoord);

  H5Fclose(file_id);
  
  return;

}

