#include "treepm_header.h"


namespace ParticleSimulator{
    namespace ParticleMesh{


static int swapShort( const short s){

  union{
    short sval;
    char cval[2];
  }u;

  u.sval = s;

  char tmp;
  tmp = u.cval[0];
  u.cval[0] = u.cval[1];
  u.cval[1] = tmp;

  return u.sval;

}


typedef struct BmpHeader{

  unsigned int  bfSize;
  unsigned short bfReserved1;
  unsigned short bfReserved2;
  unsigned int  bfOffBits;

  unsigned int  biSize;
  unsigned int  biWidth;
  unsigned int  biHeight;
  unsigned short biPlanes;
  unsigned short biBitCount;
  unsigned int  biCompression;
  unsigned int  biSizeImage;
  unsigned int  biXPixPerMeter;
  unsigned int  biYPixPerMeter;
  unsigned int  biClrUsed;
  unsigned int  biClrImporant;

}BmpHeader, *pBmpHeader;



static void outputBmp(const int width, const int height, 
		      unsigned char *color_array, const char *ofile){

  BmpHeader bmp;

  char bfType[2];
  bfType[0]          = 'B';
  bfType[1]          = 'M';
  bmp.bfSize         = width*height*3 + 54;
  bmp.bfReserved1    = 0;
  bmp.bfReserved2    = 0;
  bmp.bfOffBits      = 54;
  bmp.biSize         = 40;
  bmp.biWidth        = width;
  bmp.biHeight       = height;
  bmp.biPlanes       = 1;
  bmp.biBitCount     = 24;
  bmp.biCompression  = 0;
  bmp.biSizeImage    = 0;
  bmp.biXPixPerMeter = 0;
  bmp.biYPixPerMeter = 0;
  bmp.biClrUsed      = 0;
  bmp.biClrImporant  = 0;

#ifdef REVERSE_ENDIAN_OUTPUT
  bmp.bfSize         = swapInt(bmp.bfSize);
  bmp.bfOffBits      = swapInt(bmp.bfOffBits);
  bmp.biSize         = swapInt(bmp.biSize);
  bmp.biWidth        = swapInt(bmp.biWidth);
  bmp.biHeight       = swapInt(bmp.biHeight);
  bmp.biPlanes       = swapShort(bmp.biPlanes);
  bmp.biBitCount     = swapShort(bmp.biBitCount);
#endif

  FILE *outstream = fopen( ofile, "wb");
  fwrite( &bfType, sizeof(char), 2, outstream);
  fwrite( &bmp, sizeof(BmpHeader), 1, outstream);
  fwrite( color_array, sizeof(unsigned char), width*height*3, outstream);
  fclose( outstream);

}



static float convertLinear( const float x, const float y, const float slope,
			    float val){

  val = slope * ( val - x) + y;

  return val;

}





static float convertLog( float val){

  val = log( val + 1.0);

  return val;

}






int outputImage( const Particle *particle, const RunParam *this_run, 
		 const char *output_image){


  double range[3][2] = { {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0} };  // comoving box

  int width = IMAGEWIDTH;
  int height = IMAGEHEIGHT;
  int depth = 1;
  int npart = this_run->npart;

  double dwinv = (double)width  / (range[0][1] - range[0][0]);
  double dhinv = (double)height / (range[1][1] - range[1][0]);
  double ddinv = (double)depth  / (range[2][1] - range[2][0]);

  static float column_density[IMAGEWIDTH*IMAGEHEIGHT];
  static float column_density_sum[IMAGEWIDTH*IMAGEHEIGHT];

  int i, j;
  for( i=0; i<width*height; i++){
    column_density[i] = 0.0;
    column_density_sum[i] = 0.0;
  }

  for( i=0; i<npart; i++){
    float pos[3];
    getPos2( &particle[i], pos);
    int x = (int)(pos[0] * dwinv);
    int y = (int)(pos[1] * dhinv);
    int z = (int)(pos[2] * ddinv);
    if( x<0 || y<0 || z<0 || x>=width || y>=height || z>=depth)  continue;
#ifdef UNIFORM
    column_density[x*width+y] += this_run->uniform_mass / 2.668043e-10;
#else
    column_density[x*width+y] += particle[i].mass / 2.668043e-10;
#endif
    
  }


  // make a image on only root node 
  int nnode,inode;
  MPI_Comm_size(this_run->MPI_COMM_INTERNAL,&nnode);
  MPI_Comm_rank(this_run->MPI_COMM_INTERNAL,&inode);  
  int error_class = MPI_Reduce( column_density, column_density_sum, width*height,
				MPI_FLOAT, MPI_SUM, 0, this_run->MPI_COMM_INTERNAL);
  mpiCheck( error_class, inode);
  if( inode != 0)  return 0;

  // only inode0
  for( i=0; i<width*height; i++){
    column_density_sum[i] = convertLog( column_density_sum[i]);
    column_density_sum[i] = convertLinear( IMAGEFACA, IMAGEFACB, IMAGEFACS, 
					   column_density_sum[i]);    // refer param.h
    if( column_density_sum[i] > 255.0)  column_density_sum[i] = 255.0;
    if( column_density_sum[i] < 0.0)    column_density_sum[i] = 0.0;
  }

  unsigned char *color_array = (unsigned char *) my_malloc( sizeof(unsigned char) * width * height * 3);

  // set color map 
  int offset = 0;
  for( i=0; i<height; i++){
    for( j=0; j<width; j++){
      int index = j*height + i;
      if( COLORMAP == 2){
	float ftmp1 = -(column_density_sum[index]-255.0)*(column_density_sum[index]-255.0)/255.0 + 255.0;
	if( ftmp1 > 255.0)  ftmp1 = 255.0;
	if( ftmp1 < 0)  ftmp1 = 0.0;
	float ftmp2 = 1.3*column_density_sum[index] - 76;
	if( ftmp2 > 255.0)  ftmp2 = 255.0;
	if( ftmp2 < 0)  ftmp2 = 0.0;
	float ftmp3 = 2.0*column_density_sum[index] - 255;
	if( ftmp3 > 255.0)  ftmp3 = 255.0;
	if( ftmp3 < 0)  ftmp3 = 0.0;
	color_array[3*offset+0] = (unsigned char)ftmp3;
	color_array[3*offset+1] = (unsigned char)ftmp2;
	color_array[3*offset+2] = (unsigned char)ftmp1;
      }
      else if( COLORMAP == 1){
	float ftmp1 = 5.0*column_density_sum[index] - 1000.0;
	if( ftmp1 > 255.0)  ftmp1 = 255.0;
	if( ftmp1 < 0)  ftmp1 = 0.0;
	float ftmp2 = 1.2*column_density_sum[index] - 51;
	if( ftmp2 > 255.0)  ftmp2 = 255.0;
	if( ftmp2 < 0)  ftmp2 = 0.0;
	color_array[3*offset+0] = (unsigned char)column_density_sum[index];
	color_array[3*offset+1] = (unsigned char)ftmp2;
	color_array[3*offset+2] = (unsigned char)ftmp1;
      }
      else{
	color_array[3*offset+0] = (unsigned char)column_density_sum[index];
	color_array[3*offset+1] = (unsigned char)column_density_sum[index];
	color_array[3*offset+2] = (unsigned char)column_density_sum[index];
      }
      offset ++;
    }
  }

  outputBmp( width, height, color_array, output_image);

  outputBmp( width, height, color_array, NEWEST_IMAGEFILE);

  free(color_array);

  return 1;

}




    } // namespace ParticleMesh
}     // namespace ParticleSimulator
