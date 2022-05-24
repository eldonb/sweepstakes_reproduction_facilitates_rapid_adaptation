\pdfoutput=1
\documentclass[a4paper,12pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[lf]{Baskervaldx}
\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage[round,numbers,super]{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\setstretch{1.3}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}
\newcommand{\im}{\ensuremath{\imath} }
\newcommand{\jm}{\ensuremath{\jmath} }
\newcommand{\be}{\begin{equation}}
\newcommand{\ee}{\end{equation}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \textbf{\LARGE \@@title}

  \@@author
\end{flushright}\egroup
}
\makeatother
\title{Random sweepstakes and viability selection}
\author{ Bjarki Eldon\footnote{MfN Berlin, Germany}\footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
  through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by the Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge,   WS, and BE; acknowledge Start-up module grants through SPP 1819  with Jere Koskela and Maite Wilke-Berenguer, and  with Iulia Dahmer. \\ \today}   \orcidlink{https://orcid.org/0000-0001-9354-2391} }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}


\begin{abstract}
This C++ code generates excursions of the evolution of an advantageous type at a single locus in  a diploid
population partitioned into three genotypes, and in which two genetic types are segregating, one with viability
weight $W = \exp\left( -s(z - z_0)^2 \right)$ and the other with
viability weight one.  The population evolves according to a model of
sweepstakes reproduction and randomly occurring bottlenecks and
a simple model of viability selection. We are interested in the probability of fixation
of the advantageous type and the time to fixation conditional on
fixation, and how random sweepstakes and dominance mechanisms and
randomly occurring bottlenecks act on these statistics. 
\end{abstract}




@* {\bf Copyright}. 

Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.




@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library. 
Using a Makefile can be helpful, calling this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        c++ -Wall -Wextra -pedantic -O3 -march=native -m64 iguana.c -lm -lgsl -lgslcblas \\
        
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v --leak-check=full --show-leak-kinds=all <program call>}


@* {\bf introduction}. 
\label{SEC:intro}

We consider the evolution of an advantageous genetic type at a single
locus in a diploid population evolving according to random sweepstakes
and randomly occurring bottlenecks. 

We randomly pair individuals and each pair independently   produces  a random number of juveniles with
 probability mass function
\begin{equation}
\label{eq:1}
\prb{X_1 = k} =  \frac{ \one{2 \le k \le \Psi_N} }{ 2^{-\alpha} -   (\Psi_N + 1)^{-\alpha}  }   \left( \frac{1}{k^\alpha} - \frac{1}{(k+1)^\alpha}  \right)
\end{equation}
Write $X_1, \ldots, X_N \sim L(\alpha, \Psi_N)$ if $X_1, \ldots, X_N$ are independent and identically distributed  according to Eq~\eqref{eq:1}. Let $1 < \alpha_1 < 2 $ and $\alpha_2 > 2$ be fixed and  consider, with $\Psi_N$ fixed,  the mixture distribution
\begin{equation}
\label{eq:3}
\begin{cases}
X_1, \ldots, X_N \sim L(\alpha_1, \Psi_N) & \textrm{with probability $\varepsilon_N$,} \\
X_1, \ldots, X_N \sim L(\alpha_2, \Psi_N) & \textrm{with probability $1 - \varepsilon_N$.}
\end{cases}
\end{equation}
A scaling of $\varepsilon_N$ can be identified so that the model converges to a non-trivial limit\cite{chetwyn-diggle_beta,dahmer_coales} .

Each pair produces at least two juveniles.  Each juvenile is then
assigned two alleles, one from each parent, sampled independently and
uniformly at random. The genotype, denoted $g$,  of a given juvenile is translated into a phenotype, denoted $z$,  by a given map $f(g)$,  and  a viability weight (fitness value) $w$   is computed  
\begin{equation}
\label{eq:2}
w = \exp(-s( z_0 -  f(g) ) ) 
\end{equation}
where $z_0$ is the optimal phenotype and $s > 0$ is the fixed strength of selection.  There are obviously many ways
to configure this, and this is one of them and corresponds to the
traditional way of computing a fitness from a genotype or genotypes.
For each juvenile we sample an independent exponential with rate the
given viability weight, and the $2N$ juveniles with the smallest
exponentials survive and replace the parents.


We introduce randomly occurring bottlenecks.  In any given generation
there is a fixed probability of a bottleneck, and the size of the
bottleneck, denoted $N_b$, is fixed. When a bottleneck occurs we
randomly sample among the current individuals $N_b$ individuals to
survive the bottleneck.  The surviving individuals then produce
juveniles; if the total number of juveniles is less than the given
capacity all the juveniles survive, otherwise we assign weights and
exponentials and sample the juveniles as just described.  If all the
surviving individuals are of the optimal type we consider a fixation
has occurred since we exclude mutation, and juveniles inherit the type
of the parents as described.


Our motivation for including bottlenecks comes partly from a study on
the effect of bottlenecks on the evolution of antibiotic
resistance\cite{Mahrt2021}. The effect of selection gene genealogies
has investigated mathematically mostly in haploid populations, but recently
also in diploid populations, in particular where heterozygotes have been assumed to
be at a selective disadvantage to both
homozygotes\cite{etheridge2020genealogies}. 

@* {\bf code}.
\label{SEC:code}


In \S~\ref{SEC:struct} we collect the main containers, i.e.\ the
number of diploid individuals of each of the three possible genotypes,
the juveniles as a pair of genotype and viability weight, and the
cumulative density function for sampling the random number of
juveniles according to Eq~\eqref{eq:1} using the inverse CDF
method. We use the GSL random number generator \S~\ref{SEC:rng}, in
\S~\ref{SEC:parenttype} we sample the type of a parent, in
\S~\ref{SEC:typejuvenile} we assign a genotype to a juvenile given the
genotype of both parents, in \S~\ref{SEC:randomnumberj} we sample a
random number of juveniles using the inverse CDF method, in
\S~\ref{SEC:weight} we compute the viability weight, or the rate for
the exponential for sorting the juveniles; in \S~\ref{SEC:addjuvs} we
add juveniles produced by a parent pair with given genotypes to the
pool of juveniles; in \S~\ref{SEC:pool} we generate a pool of
juveniles for a randomly arranged pairs of parents, if there is an odd
number of individuals we assume one individual will be left out; in
\S~\ref{SEC:bottleneck} we sample number of individuals of each type
according to the hypergeometric that will survive a bottleneck; in
\S~\ref{SEC:onestep} we take one step through an excursion by checking
if a bottleneck occurs, and sampling juveniles; in
\S~\ref{SEC:trajectory} we generate one excursion; in
\S~\ref{SEC:mass} the probability mass function Eq~\eqref{eq:1} is
defined; in \S~\ref{SEC:initcdf} we compute the corresponding CDF; and
in \S~\ref{SEC:main} we generate a specified number of excursions or
experiments. 



@*1 {\bf includes}.
\label{SEC:includes}

the included libraries, and the header file containing the settings as global constants

@<includes@>=@#
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



@*1 {\bf struct M}.
\label{SEC:struct}


@<struct@>=@#
struct M {

private :

  /* diploid single locus: $0/0 = 0$; $0/1 = 1$; $1/1 = 2$; the population is partitioned into the three genotypes, the vector |population| records the number of each type indexed in the obvious way\newline */
  std::vector<double> population {} ;

/* for each juvenile we record the weight and the genotype \newline */
  std::vector< std::pair< size_t, double>> juveniles {} ;

/* the vectors containing the cumulative density functions for the random number of juveniles Eq~\eqref{eq:1} \newline */
  std::vector<double> cdf_one {} ;
  std::vector<double> cdf_two {} ;

  int numer {} ;

public :

  void setnumer (const int x) {numer = x ; }
  int getnumer(){ return numer ; }

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
  /* starting with one copy of $0/1$, the rest is homozygous $0/0$ \newline */
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
    population[0] = GLOBAL_CONST_I - 1 ;
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


@*1 {\bf the random number generator}.
\label{SEC:rng}

the GSL random number generator


@<gslrng@>=@#
gsl_rng * rngtype ; 
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


@*1 {\bf sample type of parent}. 
\label{SEC:parenttype}

Sample the type of a parent by sampling the marginal of the
multivariate hypergeometric.  There are only three possible genotypes
to sample from, and we need to sample only one copy. 

@<sample parent type@>=@#
static int sample_parent_type_mvhyper(M& self,  gsl_rng *r)
{
  int i = 0 ;
  /* the |nothers| is the total number of individuals of the type we are not sampling from \newline */
  unsigned int nothers = self.current_number_individuals() - static_cast<unsigned int>( self.number_type(i) ) ;
  unsigned int x = ( self.number_type(i) > 0 ? gsl_ran_hypergeometric( r, static_cast<unsigned>( self.number_type(i)), nothers, 1) : 0);

  if( x < 1 ){
    ++i ;
    nothers -= static_cast<unsigned int>( self.number_type(i) ) ;
    x = self.number_type(i) > 0 ? gsl_ran_hypergeometric( r, static_cast<unsigned>( self.number_type(i)), nothers, 1) : 0 ;}

/* set the sampled genotype to 2 if 0 or 1 not sampled \newline */
  i += (x < 1 ? 1 : 0);
  /* update the count of  the sampled  type in the population \newline */
  self.update_count_type( i, 1);

/* return the type of the parent \newline */
  return i ;
}


@*1 {\bf assign genotype to juvenile}.
\label{SEC:typejuvenile}


Assign genotype to a juvenile given the genotype of the two parents;
the two alleles  are sampled independently and  uniformly at random, one allele from each parent. 

@<typejuvenile@>=@#
static int assign_type_juvenile( const size_t gone, const size_t gtwo, gsl_rng *r)
{

/* |gone| and |gtwo| are the two parent genotypes \newline */
  int g {} ;
  const double u = gsl_rng_uniform(r) ;
  switch(gone){
  case 0 : {
  /* a parent genotype is $0/0$ \newline */
    g = (gtwo < 1 ? 0 : (gtwo < 2 ? (u < 0.5 ? 0 : 1) : 1) );
    break ;}
  case 1 : {
    /* a parent genotype is $0/1$ \newline */
    g = (gtwo < 1 ? (u < .5 ? 0 : 1) : (gtwo < 2 ? (u < 0.25 ? 0 : ( u < 0.75 ? 1 : 2)) : (u < 0.5 ? 1 : 2) ) ) ;
    break ; }
  case 2 : {
    /* a parent genotype is $1/1$ \newline */
    g = (gtwo < 1 ? 1 : (gtwo < 2 ? (u < .5 ? 1 : 2) : 2) ) ;
    break ; }
  default : break ; } 

  return g ;
}


@*1 {\bf sample a random number of juveniles}. 
\label{SEC:randomnumberj}


Sample a random number of juveniles; if $F$ denotes the cumulative
distribution function and $u$ a random uniform then we return $\min\left\{ j \in \mathds{N} : F(j) > u
\right\}$

@<sample random number juveniles@>=@#
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


@*1 {\bf compute viability  weight of a juvenile}.
\label{SEC:weight}

Compute a viability weight of a juvenile; the canonical form is
$\exp\left( -sf(g) \right)$ where $f$ is the trait function and $g$ the genotype.

@<weigth@>=@#
static double computeweight(const double g, gsl_rng *r)
{
  /* $z(0/1) < z(0/0) < z(1/1)$, i.e. the heterozygote is worst off \newline */
 /*  const double gg = g < 1 ? 1 : (g < 2 ? 2 : 0) ; */
 /* linear : $z(0/0) < z(0/1) < z(1/1)$; set |gg| as |2-g| \newline */
 /*  |const double gg = static_cast<double>(2 - g) ;|  */
 /*  const double gg = (g < 1 ? 2 : (g < 2 ? 1 : 0) ) ;  */
 /* recessive :  $ 0 = z(0/0) = z(0/1) < z(1/1)$, ie the beneficial type is recessive \newline */
 
 /* $r$ :  $z(0/0) < z(0/1) = z(1/1)$ ; set |delta| as |g > 0 ? 0 : 2| \newline */
 /* $j$ :  $z(0/0) = z(0/1) < z(1/1)$ ; const double gg = (g < 2 ? 2 : 0 ) ; \newline */
 /* $v$  $z(0/1) < z(0/0) < z(1/1)$  \newline */
 /* const double gg = (g < 1 ? 1 : (g < 2 ? 2 : 0) ) ; */
  const double gg = (g < 2 ? 2 : 0) ;
  /* return an exponential with rate $e^{-sf(g)}$ \newline */
  return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( gg, 2.) )) ) ;
}


@*1 {\bf adding juveniles}.
\label{SEC:addjuvs}

Add juveniles for a parent pair with given  genotypes

@<add juveniles for a given parent pair@>=@#
static void add_juveniles_for_given_parent_pair(M& self, const int gone, const int gtwo, const size_t conetwo,  gsl_rng *r)
{

/* sample the random number of juveniles \S~\ref{SEC:randomnumberj} \newline */
  const size_t numberj = sample_random_number_juveniles(self, conetwo,  r) ;
  assert( numberj > 1 ) ;
  int g {} ;
  for( size_t j = 0; j < numberj ; ++j){
  /* assign type to each juvenile using \S~\ref{SEC:typejuvenile} \newline */
    g = assign_type_juvenile( gone, gtwo, r) ;
    /* add a juvenile to the pool with weight computed in \S\ref{SEC:weight} \newline */
    self.add_juvenile( g, computeweight( static_cast<double>(g), r ) ) ; }
}


@*1 {\bf generate pool of juveniles}. 
\label{SEC:pool}

Generate a pool of juveniles

@<generate pool juveniles@>=@#
static void generate_pool_juveniles(M& self, gsl_rng *r)
{
  int gone, gtwo {} ;
  /* clear the container \newline  */
  self.removealljuveniles() ;
  /* sample distribution of number of juveniles \newline */
  const size_t conetwo = (gsl_rng_uniform(r) < GLOBAL_CONST_EPSILON ? 1 : 2) ; 
  /* |i| runs over number of pairs that can be formed from the current number of individuals \newline */
  assert( self.current_number_individuals() < GLOBAL_CONST_I + 1) ;
  const double currenti = self.current_number_individuals() ;
  for ( double i = 0 ; i < floor( currenti / 2. ) ; ++i){
  /* sample the type of the two parents \S~\ref{SEC:parenttype} \newline  */
    gone = sample_parent_type_mvhyper(self, r);
    gtwo = sample_parent_type_mvhyper(self, r);
    if( (gone > -1) && (gtwo > -1) ){
    /* should always be true \newline */
    /* add juveniles for a parent pair with given genotypes \S~\ref{SEC:addjuvs} \newline */
      add_juveniles_for_given_parent_pair( self, gone, gtwo, conetwo,  r) ;}
  }

  assert( self.totalnumberjuveniles() >= static_cast<size_t>( currenti ) );
}


@*1 {\bf sample individuals surviving a bottleneck}. 
\label{SEC:bottleneck}

Sample individuals surviving a bottleneck; the bottleneck can occur at
any time, even when the population is recovering from a previous one and has not reached  full capacity.  

@<bottleneck@>=@#
static void  bottleneck(M& self, gsl_rng *r)
{
  unsigned int n0, n1, n2 {};

  /* sample number of homozygous $0/0$  surviving a bottleneck \newline */
  n0 =  gsl_ran_hypergeometric( r, static_cast< unsigned int>(self.number_type(0) ), static_cast<unsigned int>( self.current_number_individuals() - self.number_type(0) ), GLOBAL_CONST_BOTTLENECK) ;

/* sample the number of heterozygotes  $0/1$ individuals surviving a bottleneck \newline */
  n1 =   GLOBAL_CONST_BOTTLENECK - n0 > 0 ?   gsl_ran_hypergeometric(r, static_cast<unsigned int>( self.number_type(1) ), static_cast<unsigned int>( self.current_number_individuals() -  self.number_type(0) - self.number_type(1)), GLOBAL_CONST_BOTTLENECK - n0) : 0 ;

/* set the number of homozygous $1/1$ individuals surviving a bottleneck \newline */
  n2 = GLOBAL_CONST_BOTTLENECK - n0 - n1 ;

  assert( n0 + n1 + n2 == GLOBAL_CONST_BOTTLENECK );
  /* update population with surviving copies; we are tracing the number of homozygous  $1/1$ individuals \newline */

  self.updatepopulationafterbottleneck( static_cast<double>(n0), static_cast<double>(n1), (n2 < GLOBAL_CONST_BOTTLENECK ? static_cast<double>(n2) : GLOBAL_CONST_I) ) ;
}


@*1 {\bf take one step}. 
\label{SEC:onestep}

Step through one generation in the trajectory 

@<onestep@>=@#
static void onestep(M& self, gsl_rng *r )
{
  double nth {} ;
  /* check if bottleneck \newline */
  if( gsl_rng_uniform( r) < GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
    /* bottleneck occurs; sample surviving types and update population \S~\ref{SEC:bottleneck} \newline */
    bottleneck( self, r) ; }
    if ( self.number_type(1) + self.number_type(2) > 0 ){
    /* type 1 not lost \newline */
      if( self.number_type(2) < GLOBAL_CONST_I ){
	/* not all  individuals of type 2 or homozygous $1/1$, so sample juveniles \S~\ref{SEC:pool} \newline */
	generate_pool_juveniles(self, r) ;
	if( self.totalnumberjuveniles() <= GLOBAL_CONST_II )
	  {
	    /* total number of juveniles not over capacity so all survive \newline */
	    self.update_population_all_juveniles() ;
	  }
	else{
	  /* need to partial  sort juveniles and sample according to weight \newline */
	  nth = self.nthelm() ;
	  self.sample_juveniles_according_to_weight(nth) ; }
      }
      /* mutation has fixed */
    }
    /* mutation has been lost */ 
}


@*1 {\bf generate one trajectory}.
\label{SEC:trajectory}

Generate one experiment, recording the count of the optimal type. 

@<trajectory@>=@#
static void trajectory(M& self, gsl_rng *r)
{
 /* clear containers for a new trajectory \newline */
  self.init_for_trajectory() ;

  std::vector< double > excursion_to_fixation {} ;
  /* clear the containers for the trajectory \newline */
  excursion_to_fixation.clear() ;
  excursion_to_fixation.shrink_to_fit() ;
  assert( excursion_to_fixation.size() < 1 ) ;

  
  int timi = 0 ;
  while( (self.number_type(2) + self.number_type(1) > 0 ) && (self.number_type(2) < GLOBAL_CONST_I ) ){
  assert( self.current_number_individuals() > 0 ) ;
    excursion_to_fixation.push_back( self.number_type(2) / self.current_number_individuals() ) ;
    ++ timi ;
    onestep(self, r) ; }

  size_t i = 0;
  std::fstream FILE ; 
  const std::string eskra = "trajectoryDRH" + std::to_string( self.getnumer() ) + ".txt" ;
  if( self.number_type(2) > 0 ){
    /* fixation occurs so  print excursion to file \newline */
      FILE.open(eskra, std::ios::app);
      /* print the excursion to fixation to the file \newline */
      for( const auto& y: excursion_to_fixation){ FILE << y << ' ' ;} ;
	++i ;}
        FILE << '\n' ;
      FILE.close() ;

/* print out indicator of loss or fixation and the time of the trajectory \newline */
  std::cout << (self.number_type(2) < 1 ? 0 : 1) << ' ' << timi << '\n' ;
}


@*1 {\bf the mass function}.
\label{SEC:mass}

Define the probability mass function for the random number of juveniles Eq~\eqref{eq:1}

@<mass function@>=@#
static double px(const double k, const double calpha, const double ccutoff)
{
  return ( (pow( 1./k, calpha) - pow( 1./(k + 1.), calpha) )/( pow( .5, calpha) - pow( 1./(ccutoff + 1.), calpha) ) ) ;
}



@*1 {\bf initialise CDF}.
\label{SEC:initcdf}


initialise the cumulative density function for sampling  random number of juvenile

@<initialise cdf@>=@#
static void initialise_cdf(M& self)
{
  
  for( double i = 2; i <=  GLOBAL_CONST_PSI_ONE ; ++i){
    self.add_to_cdfone( px( i, GLOBAL_CONST_ALPHA_ONE, GLOBAL_CONST_PSI_ONE) ) ;}

  for( double j = 2; j <= GLOBAL_CONST_PSI_TWO; ++j){
    self.add_to_cdftwo( px( j, GLOBAL_CONST_ALPHA_TWO, GLOBAL_CONST_PSI_TWO)) ; }

}




@*1 {\bf the main module}. 
\label{SEC:main}

@C

@<includes@>@#
@<struct@>@#
@<gslrng@>@#
@<sample parent type@>@#
@<typejuvenile@>@#
@<sample random number juveniles@>@#
@<weigth@>@#
@<add juveniles for a given parent pair@>@#
@<generate pool juveniles@>@#
@<bottleneck@>@#
@<onestep@>@#
@<trajectory@>@#
@<mass function@>@#
@<initialise cdf@>@#


int main(int argc, char *argv[])
{
  M d {} ;

  d.init_containers() ;
  initialise_cdf( d) ;
   setup_rng( static_cast<unsigned long>(atoi(argv[1]) ) );
  int z = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1;
  d.setnumer( atoi(argv[1]) ) ;
  while( --z > 0){
    trajectory(d, rngtype) ;} 
  
  d.freememory() ;

  gsl_rng_free( rngtype) ;
  
        return GSL_SUCCESS ;
}


@* {\bf examples of results}. 
\label{sec:examples}


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa}
\caption{examples of excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 0.1$, selection strength $s = 1$, bottleneck $10^2$ and probability of a bottleneck in a given generation $0.01$;  results from $10^2$ experiments with $42$ fixations and average time $57.9$ and  stdev $12.8$  }
\label{fig:ex1}
\end{figure}

\clearpage
\pagebreak
\newpage



\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufaB}
\caption{examples of excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 0.1$, selection strength $s = 1$, bottleneck $10^2$ and probability of a bottleneck in a given generation $0.1$;  results from $10^2$ experiments with $3$ fixations and average time $100.7$ and  stdev $99.6$  }
\label{fig:ex1}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufaC}
\caption{examples of excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 1$, bottleneck $10^2$ and probability of a bottleneck in a given generation $0.01$;  results from $10^2$ experiments with $31$ fixations and average time $85.4$ and  stdev $21.0$  }
\label{fig:ex1}
\end{figure}

\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufaD}
\caption{examples of excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 1$, bottleneck $10^2$ and probability of a bottleneck in a given generation $0.1$;  results from $10^2$ experiments with one fixation and average time $112$ and  stdev $--$  }
\label{fig:ex1}
\end{figure}

\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_E}
\caption{examples of excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.01$;    the dominance mechanism is linear with $z = g$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;  results from $10^2$ experiments with $46$ fixations and average time $83.8$ and  stdev $12.5$  }
\label{fig:ex1}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_F}
%%  const double gg = static_cast<double>(2 - g) ;
\caption{examples of excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.1$;  the dominance mechanism is linear with $z = g$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;   results from $10^2$ experiments with five fixations and average time $299.4$ and  stdev $187.2$  }
\label{fig:ex1}
\end{figure}



\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_G}
\caption{Examples of  excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.1$;  the dominance mechanism is $r$-shape  with $z = 2\one{g > 0}$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;   ; results from $10^2$ experiments with ten fixations and average time to fixation  $1600.2$ and stdev  $631.3$ }
\label{fig:excursprufaG}
\end{figure}


\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_H}
\caption{Examples of  excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 0.1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.1$;  the dominance mechanism is $r$-shape  with $z = 2\one{g > 0}$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;   ; results from $10^2$ experiments with one fixation }
\label{fig:exampleH}
\end{figure}

\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_I}
\caption{Examples of  excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 0.1$, selection strength $s = 0.1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.01$;  the dominance mechanism is $r$-shape  with $z = 2\one{g > 0}$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;   ; results from $10^2$ experiments with ten fixations and average time to fixation $818.8$ and stdev $191.6$ }
\label{fig:exampleI}
\end{figure}

\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_J}
\caption{Examples of  excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 0.1$, bottleneck $10^2$ and probability of a bottleneck in a given generation $0.01$;  the dominance mechanism is $r$-shape  with $z = 2\one{g > 0}$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;   ; results from $10^2$ experiments  ; fourteen fixations with 801.0000 380.229 }
\label{fig:exampleJ}
\end{figure}

\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_K}
\caption{Examples of  excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 1.05$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 0.1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.1$;  the dominance mechanism is $r$-shape  with $z = 2\one{g > 0}$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$;  eleven fixations out of $10^3$ experiments with average time to fixation $4692.2$ and stdev $2403.0$}
\label{fig:exampleK}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_L}
\caption{Examples of  excursions to fixation of the homozygous type $1/1$ at a single locus in a diploid population evolving according to random sweepstakes Eq~\eqref{eq:1} and randomly occurring bottlenecks with  fixed capacity $2N=10^6$, cutoff $\Psi_N = 2N$, $\alpha_1 = 0.75$, $\alpha_2 = 3$,  $\varepsilon_N = 1/N$, selection strength $s = 0.1$, bottleneck $10^4$ and probability of a bottleneck in a given generation $0.1$;  the dominance mechanism is linear  with $z =    g$ for $g \in \{0,1,2\}$,  optimal type  $1/1$ or $z_0 = 2$; five fixations out of $10^3$ experiments with  average time to fixation $1397.2$ and stdev  $508.9$  }
\label{fig:exampleL}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excurs_prufa_M}
\caption{ex ; 170.0000 410.7588  49.475 }
\label{fig:exampleM}
\end{figure}


\clearpage
\pagebreak
\newpage

@* {\bf conclusion}. 
\label{sec:concl}


\cite{etheridge2020genealogies} consider a model of selection in a
diploid population in which the homozygous $aa$ individuals are at a
selective disadvantage to $AA$ individuals and heterozygous $aA$
individuals are at a selective disadvantage to both homozygotes. 

@* {\bf references}. 
\label{sec:refs}

\bibliographystyle{plain}
\bibliography{/home/bjarki/verk/master_bibfile/refs.bib}

@
\end{document}
