#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <cstddef>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <chrono>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "diploid_excursions_random_bottlenecks.hpp"

struct M {

private :

  /* diploid single locus: 0/0 = 0; 0/1 = 1; 1/1 = 2 */
  std::vector<double> population {} ;

  std::vector< std::pair< size_t, double>> juveniles {} ;

  std::vector<double> cdf_one {} ;
  std::vector<double> cdf_two {} ;

public :

  void removealljuveniles()
  { juveniles.clear() ;
    juveniles.shrink_to_fit();
    assert( juveniles.size() < 1 );
  }
    
  size_t totalnumberjuveniles()
  { return juveniles.size() ; }
  
  void add_juvenile( const size_t g, const double weight)
  { juveniles.push_back( std::make_pair(g, weight) ) ; }
  
  double number_type( const size_t c_type )
  { return population[c_type]; }

  double current_number_individuals()
  { return (population[0] + population[1] + population[2]) ; }

  void updatepopulationafterbottleneck( const double n0, const double n1, const double n2)
  {
    population[0] = n0 ;
    population[1] = n1 ;
    population[2] = n2 ;
    
  }

  void update_count_type( const size_t c_type, const size_t add_subtract )
  { population[c_type ] +=  (add_subtract < 1 ? 1. : -1.) ; }

  size_t type_fewest_copy()
  { return std::distance( population.begin(), std::min_element(population.begin(), population.end()) );  }

  size_t type_most_copies()
  { return std::distance( population.begin(), std::max_element(population.begin(), population.end()) );  }

  void update_population_all_juveniles()
  {
    std::fill( population.begin(), population.end(), 0. );
    for( const auto &j : juveniles){
      population[ std::get<0>( j ) ] += 1;}
  }

  double get_cdfone(const size_t i)
  { return cdf_one[i] ; }

  double get_cdftwo(const size_t i)
  { return cdf_two[i];}

  void add_to_cdfone( const double c_x)
  { cdf_one.push_back( cdf_one.back() +  c_x ) ; }

  void add_to_cdftwo( const double c_x)
  { cdf_two.push_back( cdf_two.back() + c_x ) ; }

  void init_containers()
  {
    population.reserve(3);
    population[0] = GLOBAL_CONST_I - 1 ;
    population[1] = 1 ;
    population[2] = 0 ;
    cdf_one.reserve( GLOBAL_CONST_CUTOFF_ONE + 2) ;
    cdf_two.reserve( GLOBAL_CONST_CUTOFF_TWO + 2) ;
    assert( population[0] + population[1] + population[2] == GLOBAL_CONST_I);
    cdf_one.push_back(0.);
    cdf_one.push_back(0.);
    cdf_two.push_back(0.);
    cdf_two.push_back(0.);
    assert( cdf_one.size() == 2);
    assert( cdf_two.size() == 2);
  }

  void init_for_trajectory()
  {
    population[0] = GLOBAL_CONST_I -1 ;
    population[1] = 1 ;
    population[2] = 0 ;
    assert( population[0] + population[1] + population[2] == GLOBAL_CONST_I ) ;
  }

 static bool comp( std::pair<size_t, double> a, std::pair<size_t, double> b) 
  { 
    return ( std::get<1>(a) < std::get<1>(b) ); 
  }

  double nthelm( )
{
  std::nth_element( juveniles.begin(), juveniles.begin() + (GLOBAL_CONST_II - 1), juveniles.end(), comp);

  return( std::get<1>(juveniles[ GLOBAL_CONST_II - 1]) );
}

  void sample_juveniles_according_to_weight(const double c_nth)
  {
    population[0] = 0 ;
    population[1] = 0 ;
    population[2] = 0 ;
    size_t j = 0 ;
    while( j < GLOBAL_CONST_II ){
      assert( j < GLOBAL_CONST_II) ;
      population[ std::get<0>(juveniles[j]) ] += std::get<1>(juveniles[j]) <= c_nth ? 1 : 0 ;
      ++j ;
    }
    assert( population[0] + population[2] + population[1] == GLOBAL_CONST_I  );
  }
  
  void freememory()
  {
    population.clear();
    std::vector<double>().swap(population);
    juveniles.clear() ;
    std::vector< std::pair<size_t, double> >().swap( juveniles) ;
    cdf_one.clear() ;
    std::vector<double>().swap( cdf_one) ;
    cdf_two.clear() ;
    std::vector<double>().swap( cdf_two) ;
  }
} ;

gsl_rng * rngtype ; 
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


static int sample_parent_type(M& self, gsl_rng *r)
{

  size_t m_i = self.type_fewest_copy() ;
  size_t m_a = self.type_most_copies() ;
  double u = gsl_rng_uniform(r) ;
  int g {} ;
  if( self.current_number_individuals() > 1 ){
  switch(m_i){
  case 0 : {
    if( m_a > 1 ){
      /* order is 0, 1, 2 */
      g =  (u < self.number_type(0)/( self.current_number_individuals() ) ? 0 : ( u < (self.number_type(0) + self.number_type(1))/self.current_number_individuals() ? 1 : 2)) ;
      /* update count of types */
      self.update_count_type( g, 1);}
    else{
      /* order is 0, 2, 1 */
      g =  (u < self.number_type(0)/( self.current_number_individuals() ) ? 0 : ( u < (self.number_type(0) + self.number_type(2))/self.current_number_individuals() ? 2 : 1)) ;
       /* update count of types */
      self.update_count_type( g, 1);}
    break ; }
  case 1 : {
    if( m_a > 0 ){
      /* order is 1, 0, 2 */
      g =  (u < self.number_type(1)/( self.current_number_individuals() ) ? 1 : ( u < (self.number_type(1) + self.number_type(0))/self.current_number_individuals() ? 0 : 2)) ;
        /* update count of types */
      self.update_count_type( g, 1);}
    else{
      /* order is 1, 2, 0 */
      g =  (u < self.number_type(1)/( self.current_number_individuals() ) ? 1 : ( u < (self.number_type(1) + self.number_type(2))/self.current_number_individuals() ? 2 : 0)) ;
        /* update count of types */
      self.update_count_type( g, 1);}
    break ;  }
  case 2 : {
    if (m_a > 0){
      /* order is 2, 0, 1 */
      g =  (u < self.number_type(2)/( self.current_number_individuals() ) ? 2 : ( u < (self.number_type(2) + self.number_type(0))/self.current_number_individuals() ? 0 : 1)) ;
        /* update count of types */
      self.update_count_type( g, 1);}
    else{
      /* order is 2, 1, 0 */
      g =  (u < self.number_type(2)/( self.current_number_individuals() ) ? 2 : ( u < (self.number_type(2) + self.number_type(1))/self.current_number_individuals() ? 1 : 0)) ;
        /* update count of types */
      self.update_count_type( g, 1);}
    break ; }
  default : break ; }
  }
  else{
    /* not more than one individual to choose from */
    g = (self.current_number_individuals() < 1 ? -1 : self.type_most_copies()) ;
    self.update_count_type( g, 1);
  }
  
  return g ;
}


static int sample_parent_type_mvhyper(M& self,  gsl_rng *r)
{
  int i = 0 ;
  unsigned int nothers = self.current_number_individuals() - static_cast<unsigned int>( self.number_type(i) ) ;
  unsigned int x = gsl_ran_hypergeometric( r,  self.number_type(i), nothers, 1);

  if( x < 1 ){
    ++i ;
    nothers -= static_cast<unsigned int>( self.number_type(i) ) ;
    x = gsl_ran_hypergeometric( r,  self.number_type(i), nothers, 1);}

  i += (x < 1 ? 1 : 0);
  self.update_count_type( i, 1);

  return i ;
}



/* assign genotype to juvenile given parent genotype */
static int assign_type_juvenile( const size_t gone, const size_t gtwo, gsl_rng *r)
{

  int g {} ;
  const double u = gsl_rng_uniform(r) ;
  switch(gone){
  case 0 : {
    g = (gtwo < 1 ? 0 : (gtwo < 2 ? (u < 0.5 ? 0 : 1) : 1) );
    break ;}
  case 1 : {
    g = (gtwo < 1 ? (u < .5 ? 0 : 1) : (gtwo < 2 ? (u < 0.25 ? 0 : ( u < 0.75 ? 1 : 2)) : (u < 0.5 ? 1 : 2) ) ) ;
    break ; }
  case 2 : {
    g = (gtwo < 1 ? 1 : (gtwo < 2 ? (u < .5 ? 1 : 2) : 2) ) ;
    break ; }
  default : break ; } 

  return g ;
}

static size_t sample_random_number_juveniles(M& self, const size_t c_twoone,  gsl_rng *r)
{
  const double u = gsl_rng_uniform(r);
  size_t j = 2 ;
  if( c_twoone < 2 ){
    while( u > self.get_cdfone(j) ){ ++j ;}}
  else{
    while( u > self.get_cdftwo(j) ){ ++j ;}}

  assert( j > 1) ;
  return j ;
}

static double computeweight(const double g, gsl_rng *r)
{
    /*   z(0/0) < z(0/1) = z(1/1) : ie the fitness jumps to top fitness at 0/1
   ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow(2. - ( g > 0 ? 2. : 0.), 2.) )) ) ; */
  /* :  z(0/0) < z(0/1) < z(1/1), ie the fitness is linearly increasing  */
  /*   ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow(2. - ( g ), 2.) )) ) ; */
  /* now j : z(0/0) = z(0/1) < z(1/1) */
  /* return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( g < 2 ? 2 : 0, 2.) )) ) ; */
  /* z(0/1) < z(0/0) < z(1/1), i.e. the heterozygote is worst off */
  const double gg = g < 2 ? 2 : 0 ;
  return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( gg, 2.) )) ) ;
  /*
  const double gg = (g < 1 ? 1 : (g < 2 ? 2 : 0) ) ;
  return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( gg, 2.) )) ) ;
  */

}


static void add_juveniles_for_given_parent_pair(M& self, const int gone, const int gtwo, const size_t conetwo,  gsl_rng *r)
{

  const size_t numberj = sample_random_number_juveniles(self, conetwo,  r) ;
  assert( numberj > 1 ) ;
  int g {} ;
  for( size_t j = 0; j < numberj ; ++j){
    g = assign_type_juvenile( gone, gtwo, r) ;
    self.add_juvenile( g, computeweight( static_cast<double>(g), r ) ) ; }
}


static void generate_pool_juveniles(M& self, gsl_rng *r)
{
  int gone, gtwo {} ;
  self.removealljuveniles() ;
  /* sample distribution of number of juveniles */
  const size_t conetwo = (gsl_rng_uniform(r) < GLOBAL_CONST_EPSILON ? 1 : 2) ; 
  /* i runs over number of pairs that can be formed from the current number of individuals */
  assert( self.current_number_individuals() < GLOBAL_CONST_I + 1) ;
  const double currenti = self.current_number_individuals() ;  
  for ( double i = 0 ; i < floor( currenti / 2. ) ; ++i){
    gone = sample_parent_type(self, r);
    gtwo = sample_parent_type(self, r);
    if( (gone > -1) && (gtwo > -1) ){
      add_juveniles_for_given_parent_pair( self, gone, gtwo, conetwo,  r) ;}
  }
  /*  printf( "j  %g   %lu\n",  currenti, self.totalnumberjuveniles() ) ; */
  assert( self.totalnumberjuveniles() >= static_cast<size_t>( currenti ) );
 
}



/* sample number of each type 0/1/2 surviving a bottleneck */
static void  bottleneck(M& self, gsl_rng *r)
{
  /* sample number of type 0 surviving a bottleneck */
  unsigned int n0, n1, n2 {};
  
  n0 =  gsl_ran_hypergeometric( r, static_cast< unsigned int>(self.number_type(0) ), static_cast<unsigned int>( self.current_number_individuals() - self.number_type(0) ), GLOBAL_CONST_BOTTLENECK) ;

  n1 =   GLOBAL_CONST_BOTTLENECK - n0 > 0 ?   gsl_ran_hypergeometric(r, static_cast<unsigned int>( self.number_type(1) ), static_cast<unsigned int>( self.current_number_individuals() -  self.number_type(0) - self.number_type(1)), GLOBAL_CONST_BOTTLENECK - n0) : 0 ;

  n2 = GLOBAL_CONST_BOTTLENECK - n0 - n1 ;

  assert( n0 + n1 + n2 == GLOBAL_CONST_BOTTLENECK );
  /* update population with surviving copies */
  //printf("bf  %g %g %g\n", self.number_type(0), self.number_type(1),self.number_type(2) ) ;
  self.updatepopulationafterbottleneck( static_cast<double>(n0), static_cast<double>(n1), (n2 < GLOBAL_CONST_BOTTLENECK ? static_cast<double>(n2) : GLOBAL_CONST_I) ) ;
  // printf("ba  %u %u %u\n", n0, n1, n2) ;
}


static void onestep(M& self, gsl_rng *r )
{
  double nth {} ;
  /* check if bottleneck */
  if( gsl_rng_uniform( r) < GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
    /* bottleneck occurs ; sample surviving types and update population */
    bottleneck( self, r) ; }
    if ( self.number_type(1) + self.number_type(2) > 0 ){
      if( self.number_type(2) < GLOBAL_CONST_I ){
	/* not all  individuals of type 2, so sample juveniles */
	generate_pool_juveniles(self, r) ;
	if( self.totalnumberjuveniles() <= GLOBAL_CONST_II )
	  {
	    /* total number of juveniles not over capacity so all survive */
	    self.update_population_all_juveniles() ;
	  }
	else{
	  /* need to sort juveniles and sample according to weight */
	  nth = self.nthelm() ;
	  self.sample_juveniles_according_to_weight(nth) ; }
      }
      /* mutation has fixed */
    }
    /* mutation has been lost */ 
}


static void trajectory(M& self, gsl_rng *r)
{
  self.init_for_trajectory() ;

  std::vector< double > excursion_to_fixation {} ;
  std::vector< double > v_current_number_individuals {} ;
  excursion_to_fixation.clear() ;
  excursion_to_fixation.shrink_to_fit() ;
  v_current_number_individuals.clear() ;
  int timi = 0 ;
  while( (self.number_type(2) + self.number_type(1) > 0 ) && (self.number_type(2) < GLOBAL_CONST_I ) ){
    v_current_number_individuals.push_back( self.current_number_individuals() ) ;
    excursion_to_fixation.push_back( self.number_type(2)  ) ;
    ++ timi ;
    onestep(self, r) ; }

  size_t i = 0;
  if( self.number_type(2) > 0 ){
    /* fixation occurs so  print excursion to file */
      FILE *fptr= fopen("acase_excursions_diploid_random_bottlenecks_zj_s01.txt", "a");
      for( const auto& y: excursion_to_fixation){ fprintf(fptr, "%g ", y/v_current_number_individuals[i] ) ;
	++i ;}
      fprintf(fptr, "\n");
      fclose( fptr);
  }
  
  std::cout << (self.number_type(2) < 1 ? 0 : 1) << ' ' << timi << '\n' ;
}

static double px(const double k, const double calpha, const double ccutoff)
{
  return ( (pow( 1./k, calpha) - pow( 1./(k + 1.), calpha) )/( pow( .5, calpha) - pow( 1./(ccutoff + 1.), calpha) ) ) ;
}


static void initialise_cdf(M& self)
{
  
  for( double i = 2; i <=  GLOBAL_CONST_PSI_ONE ; ++i){
    self.add_to_cdfone( px( i, GLOBAL_CONST_ALPHA_ONE, GLOBAL_CONST_PSI_ONE) ) ;}

  for( double j = 2; j <= GLOBAL_CONST_PSI_TWO; ++j){
    self.add_to_cdftwo( px( j, GLOBAL_CONST_ALPHA_TWO, GLOBAL_CONST_PSI_TWO)) ; }

}


static void testing (M& self )
{

  std::vector< std::pair< size_t, double> > juvs {} ;

  juvs.push_back( std::make_pair< size_t, double>(0, .1) ) ;
  juvs.push_back( std::make_pair< size_t, double>(1, .2) ) ;
  juvs.push_back( std::make_pair< size_t, double>(3, .9) ) ;

  assert(  self.comp( juvs[0], juvs[1] ) ) ;
  assert(  self.comp( juvs[1], juvs[2] ) ) ;
  assert(  self.comp( juvs[0], juvs[2] ) ) ;

  std::cout << std::get<0>( juvs[0] ) << ' ' <<  std::get<0>( juvs[1] ) << ' ' <<  std::get<0>( juvs[2] ) <<  '\n' ;

  std::cout << self.get_cdfone(0) << ' ' << self.get_cdfone( 1) << ' ' << self.get_cdfone(2) << ' ' << self.get_cdfone(3) << ' ' <<   self.get_cdfone(GLOBAL_CONST_CUTOFF_ONE) << '\n' ;
  std::cout << self.get_cdftwo(0) << ' ' << self.get_cdftwo( 1) << ' ' << self.get_cdftwo(2) << ' ' <<  self.get_cdftwo(3) << ' ' <<  self.get_cdftwo(GLOBAL_CONST_CUTOFF_TWO) << '\n' ;
}




int main(int argc, char *argv[])
{

  /*  printf("I %g\n", GLOBAL_CONST_I) ; */
  
  M d {} ;

  d.init_containers() ;
  initialise_cdf( d) ;

   testing(d) ;
  
  setup_rng( static_cast< unsigned long>(atoi(argv[1]) ) );
  int z = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1;
  while( --z > 0){
    trajectory(d, rngtype) ;} 
  
  d.freememory() ;

  gsl_rng_free( rngtype) ;
  return 0;
}
