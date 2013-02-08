extern "C"{
#include "predicates.h"
}

int test( double eps_x, double eps_y, int exact )
{
  
  double mpa[2], mpb[2], mpc[2], mpd[2];
  mpa[0] = 1.0;
  mpa[1] = 0.0;
  mpc[0] = 1.0;
  mpc[1] = 1.0;
  mpb[0] = 0.0;
  mpb[1] = 1.0;
  mpd[0] = eps_x;
  mpd[1] = eps_y;
  
  double orientation;
  double inCircle;
  int result;

  //--------------------------------------------------------------------------------------------
  // Test Is in Circle
  if( exact )
    {
    orientation = orient2d( mpa, mpb, mpc );
    inCircle = incircle( mpa, mpb, mpc, mpd ) * orientation;
    }
  else
    {
    orientation = orient2dfast( mpa, mpb, mpc );
    inCircle = incirclefast( mpa, mpb, mpc, mpd) * orientation;
    }
  
  if( inCircle < 0 )
    {
    result = 1;
    }
  else
    {
    result = 0;
    }
  
  return result;
}

int predicates_test( int , char** )
{

  printf( "Epsilon_x: %f \n", -0.2);
  printf( "Epsilon_y: %f \n", -0.2);
  printf( "Exact Testing? %d \n" , 0);
  int test1 = test( -0.2, -0.2, 0 );
  printf("Result? %d\n", test1);
  printf("The result should fail (1)\n");

  printf( "Epsilon_x: %f \n", 0.2);
  printf( "Epsilon_y: %f \n", 0.2);
  printf( "Exact Testing? %d \n" , 0);
  int test2 = test( 0.2, 0.2, 0 );
  printf("Result? %d\n", test2);
  printf("The result should pass (0)\n");
  
  printf( "Epsilon_x: %f \n", 0.0);
  printf( "Epsilon_y: %f \n", -0.00000000000000002);
  printf( "Exact Testing? %d \n" , 0);
  int test3 = test( 0.0, -0.00000000000000002, 0 );
  printf("Result? %d\n", test3);
  printf("The result should pass (0)\n");
  
  printf( "Epsilon_x: %f \n", 0.0);
  printf( "Epsilon_y: %f \n", -0.00000000000000002);
  printf( "Exact Testing? %d \n" , 1);
  int test4 = test( 0.0, -0.00000000000000002, 1 );
  printf("Result? %d\n", test4);
  printf("The result should fail (1)\n");
  
  return 0;
}
