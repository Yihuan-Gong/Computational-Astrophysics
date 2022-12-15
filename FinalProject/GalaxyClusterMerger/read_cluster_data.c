#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#include "constants.h"
#include "mangle_names.h"

void FTOC(read_num_points)(char filename[MAX_STRING_LENGTH+1], int *num_points)
{

  hid_t   file_id, dataset, dataspace;
  herr_t  status;
  hsize_t dims[1], maxdims[1];

  int rank;
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

  dataset   = H5Dopen(file_id, "/fields/radius");
  dataspace = H5Dget_space(dataset);
  rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  *num_points = (int)dims[0];

  H5Fclose(file_id);

  return;

}

void FTOC(read_particle_number)(char filename[MAX_STRING_LENGTH+1], int *num_particles)
{

  hid_t   file_id, dataset, dataspace;
  herr_t  status;
  hsize_t dims[1], maxdims[1];

  int rank;
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

  dataset   = H5Dopen(file_id, "/dm/particle_mass");
  dataspace = H5Dget_space(dataset);
  rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  *num_particles = (int)dims[0];

  H5Fclose(file_id);

  return;

}

void FTOC(read_profile)(char filename[MAX_STRING_LENGTH+1], 
			char fieldname[MAX_STRING_LENGTH+1],
			double field[])
{

  hid_t   file_id, dataset;
  herr_t  status;

  int len;
  char* string_index;
  char local_filename[MAX_STRING_LENGTH+1];
  char local_fieldname[MAX_STRING_LENGTH+1];

  len = 0;
  string_index = filename;

  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';

  len = 0;
  string_index = fieldname;

  while (*string_index != ' ') {
    local_fieldname[len] = fieldname[len];
    len++;
    string_index++;
  }

  *(local_fieldname+len) = '\0';

  file_id = H5Fopen(local_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset   = H5Dopen(file_id, local_fieldname);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                   H5S_ALL, H5P_DEFAULT, field);
  H5Dclose(dataset);

  H5Fclose(file_id);

  return;

}

void FTOC(read_profiles)(char filename[MAX_STRING_LENGTH+1], int *hydro_limit,
			 double r[], double pden[], double dens[], double pres[],
			 double gpot[], double grav[], double metl[], double bmag[])
{

  hid_t   file_id, dataset;
  herr_t  status;

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

  dataset   = H5Dopen(file_id, "/fields/radius");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, r);
  H5Dclose(dataset);

  /*
  dataset   = H5Dopen(file_id, "/fields/dark_matter_density");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, pden);
  H5Dclose(dataset);
  */

  dataset   = H5Dopen(file_id, "/fields/density");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, dens);
  H5Dclose(dataset);

  dataset   = H5Dopen(file_id, "/fields/pressure");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, pres);
  H5Dclose(dataset);

  dataset   = H5Dopen(file_id, "/fields/gravitational_potential");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, gpot);
  H5Dclose(dataset);

  dataset   = H5Dopen(file_id, "/fields/gravitational_field");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, grav);
  H5Dclose(dataset);

  /*
  dataset   = H5Dopen(file_id, "/fields/metallicity");
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		   H5S_ALL, H5P_DEFAULT, metl);
  H5Dclose(dataset);
  */

  if (!(*hydro_limit)) {
    dataset   = H5Dopen(file_id, "/fields/magnetic_field_strength");
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
		     H5S_ALL, H5P_DEFAULT, bmag);
    H5Dclose(dataset);
  }
      
  H5Fclose(file_id);

  return;

}

void FTOC(read_particle_mass)(char filename[MAX_STRING_LENGTH+1], double *pmass)
{
  hid_t   file_id, dataset, dataspace, memspace;
  herr_t  status;

  int len = 0;
  char* string_index;
  char local_filename[MAX_STRING_LENGTH+1];

  hsize_t start[1], stride[1], count[1], dims[1];
  int rank = 1;

  start[0] = 0;
  stride[0] = 1;
  count[0] = 1;
  dims[0] = 1;

  string_index = filename;

  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';

  file_id = H5Fopen(local_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset   = H5Dopen(file_id, "/dm/particle_mass");
  dataspace = H5Dget_space(dataset);
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, pmass);

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Fclose(file_id);

}

void FTOC(read_particles)(char filename[MAX_STRING_LENGTH+1], 
			  double xpos[], double ypos[], double zpos[], 
			  double xvel[], double yvel[], double zvel[])
{

  hid_t   file_id, dataset, dataspace, memspace;
  herr_t  status;

  hsize_t start[2], stride[2], count[2], dims[2], maxdims[2];

  int len = 0;
  char* string_index;
  char local_filename[MAX_STRING_LENGTH+1];

  int rank;

  string_index = filename;
  
  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';

  stride[0] = 1;
  stride[1] = 1;

  start[0] = 0;

  file_id = H5Fopen(local_filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset   = H5Dopen(file_id, "/dm/particle_position");

  dataspace = H5Dget_space(dataset);

  rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  count[0] = dims[0];
  count[1] = 1;

  dims[1] = 1;
 
  start[1] = 0;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, xpos);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 1;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, ypos);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 2;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, zpos);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);

  H5Dclose(dataset);

  dataset   = H5Dopen(file_id, "/dm/particle_velocity");

  dataspace = H5Dget_space(dataset);

  start[1] = 0;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, xvel);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 1;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, yvel);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  dataspace = H5Dget_space(dataset);

  start[1] = 2;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                               stride, count, NULL);
  memspace = H5Screate_simple(rank, dims, NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                   H5P_DEFAULT, zvel);

  H5Sclose(memspace);
  H5Sclose(dataspace);

  H5Dclose(dataset);

  H5Fclose(file_id);

  return;

}
