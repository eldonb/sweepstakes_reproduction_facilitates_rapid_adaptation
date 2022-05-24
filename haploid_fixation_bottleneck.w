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
\setstretch{1.5}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
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
\title{Fixing with a fat tail}
\author{Bjarki Eldon\footnote{MfN Berlin, Germany} \footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
  through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by the Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge,   WS, and BE. BE also acknowledges Start-up module grants through SPP 1819  with Jere Koskela and Maite Wilke-Berenguer, and  with Iulia Dahmer. \\ \today} \orcidlink{https://orcid.org/0000-0001-9354-2391} }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}


\begin{abstract}
 This code generates excursions of the evolution of a haploid
      population partitioned into two genetic types, with viability
      weight determined by $W = e^{-s(g - g_0)^2}$, where $g$ is the
      genetic type of a given individual, and $g_0$ is the optimal
      type, and $s > 0$ is the strength of selection.  The population
      evolves according to a model of random sweepstakes and viability
      selection and randomly occurring bottlenecks.  We estimate the probability of fixation of the type
      conferring advantage, and the expected time to fixation
      conditional on fixation of the advantageous type. 
\end{abstract}

\tableofcontents


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
.c} file), one needs the GNU Scientific Library, and the C++ boost library.  
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
\label{intro}


We consider a haploid population of fixed  size $N$. Let $X^N, X_1^N, \ldots, X_N^N$ be i.i.d.\ discrete random variables taking values in $\{1, \ldots, \Psi_N\}$; the $X_1^N, \ldots, X_N^N$ denote the random number of juveniles independently   produced  in a given generation according to
\be
\label{PXN}
   \prb{X^N = k} = \frac{ (\Psi_N +1)^\alpha }{ (\Psi_N + 1)^\alpha -
1 } \left( \frac{1}{k^\alpha} - \frac{1}{(k+1)^{\alpha}} \right),
\quad 1 \leq k \leq \Psi_N.  \ee   The mass in Eq \eqref{PXN} is normalised so that $\prb{ 1 \leq
X^N \leq \Psi_N} =1 $, and $\prb{X^N = k} \ge \prb{X^N = k+1}$. Given
a pool of at least $N$ juveniles, we sample $N$ juveniles for the next
generation.  Leaving out an atom at zero gives $X_1^N + \cdots + X_N^N
\ge N$ almost surely, guaranteeing that we always have at least $N$
juveniles to choose from in each generation. 


Write $X_{1} \sim L(\alpha,\Psi_{N})$ if $X_{1}$ is distributed
according to Eq~\eqref{PXN} for given values of $\alpha$ and
$\Psi_{N}$. Let $1 < \alpha_{1} < 2$ and $\alpha_{2} > 2$ be fixed and
consider the mixture distribution\cite{dahmer_coales}
\begin{equation}
\label{eq:3}
X_{1}, \ldots , X_{N} \sim
\begin{cases}
L(\alpha_{1}, \Psi_{N}) & \text{with probability $\varepsilon_{N}$,} \\
L(\alpha_{2}, \Psi_{N}) & \text{with probability $1 - \varepsilon_{N}$.} \\
\end{cases}
\end{equation}
Similarly, by identifying the appropriate scaling of $\varepsilon_{N}$
one can keep
$\alpha$ fixed and varied $\Psi_{N}$\cite{chetwyn-diggle_beta}.  

Each juvenile inherits the type of its' parent, and is assigned a
viability weight $z$ accoring to the type; the wild type is assigned
the weight $z=e^{-s}$ for some fixed $s > 0$, and  advantageous type
the weight one.  For each juvenile we sample an exponential with rate
the given viability weight, and the $N$ juveniles with the smallest
exponential replace the parents.  In any given generation a bottleneck
of a fixed size $N_{b}$
occurs with a fixed probability. If a bottleneck occurs we sample
$N_{b}$ individuals independently and uniformly at random without
replacement. The surviving individuals then produce juveniles, and if
the total number of juveniles is less than  the capacity $N$ all the
juveniles survive, otherwise we assign weights and sample $N$
juveniles according to the weights.

Let $Y_{t} \equiv \{Y_{t} : t \ge 0 \}$ denote the frequency of the
type conferring selective advantage, and write  $T_{k}(y) :=  \min\{ t
\ge 0 :  Y_{t} = k, Y_{0} = y \}$. We are interested in the
quantitites
\begin{equation}
\label{eq:4}
\begin{split}
p_{N} & :=   \prb{ T_{N}(1) < T_{0}(1) } \\
\tau_{N} & := \EE{ T_{N}(1) :  T_{N}(1) < T_{0}(1) } \\
\end{split}
\end{equation}

@* {\bf Code}. 
\label{SEC:code}

We collect the key containers and constants  into a struct
\S~\ref{SEC:structM},  we use the GSL random  number generator
\S~\ref{SEC:rng}, in \S~\ref{SEC:cdf} we compute the cumulative
density function for sampling random numbers of juveniles according to
the inverse CDF method,  in \S~\ref{SEC:random_number_juveniles} we
sample a random number of juveniles, in \S~\ref{SEC:comp} we define a
comparison function for sorting the exponentials in \S~\ref{SEC:nth},
in \S~\ref{SEC:samplepool} we sample a pool of juveniles and assign
weight to them in \S~\ref{SEC:assignweight},  in
\S~\ref{sec:survivehypergeom} we sample the number of individuals  of the advantageous
type surviving  a bottleneck,  in \S~\ref{sec:surviveweight} we count
the number of advantageous type surviving selection according to their
weight,  in \S~\ref{sec:onestep} we step through one generation by
checking if a bottleneck occurs and then produce juveniles if neither
fixation nor loss of the advantageous type occurs,  in
\S~\ref{sec:trajectory} we generate one excursion until fixation or
loss of the advantageous type starting with one copy of the
advantageous type, the main module \S~\ref{SEC:main} generates a given
number of trajectories, \S~\ref{sec:examples} holds examples of
trajectories to fixation of the advantageous type. 

@*1 {Includes}. 
\label{SEC:includes}

The included libraries.

@<includes@>=@#
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
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


@*1 {\bf the data struct}.
\label{SEC:structM}

The data structure collecting the main containers and the  constants;
we start with the type conferring advantage in one copy, unless
otherwise stated. 

@<structM@>=@#

struct M {
  const size_t N = 1000000;
  /* keep |alpha_one| less than or equal to   |alpha_two| \newline */
  const double alpha_one = 0.75 ;
  const double alpha_two = 3.0;
  /* {\it psi} is the cutoff $\Psi_N$ \newline */
  const size_t cutoff = N;
  const  size_t bottlenecksize = 10000;
  const double pbottle = -1. ;
  /* psi is psione and psione / N to 0; need $pN = (psione / N)^{2 - alpha}$  */
  const double psi =  static_cast<double>(cutoff); 
  const double psitwo =  static_cast<double>( cutoff) ;
  /* |c_strength_selection| the strength of selection $s > 0$ \newline */
  const double c_strength_selection = 0.5 ;
  /*  for  cutoffs need pN as epsilonN and  $pN =  (psione / N)^{ 2 -
  alpha }$ \newline   */
  /* |epsilon_N| is the probability of producing juveniles with |alpha_one| or  |psitwo|; consider as $\varepsilon_N = c/N$ where |c_c| is $c > 0$ a constant \newline */
    /*  for  cutoffs need pN as epsilonN and  $pN =  (psione / N)^{ 2 - alpha }$  */
  const double epsilon_N =  0.1 ;
    /*  |1./static_cast<double>(N)| ; \newline */
    
  /* probability of  psi is one minus epsilon; probability of psi2 is  epsilon */
  const int number_experiments = 250;
  int numer {} ;
  const std::string c_excursion_skra = "excursionskraparallelnobottle" ;
  size_t SN = 0;
  size_t SNW = 0;
  /* initial number of copies of advantageous type \newline */
  size_t y = 1;
  size_t Nprime {} ;
  double Nth {} ;
  /* CDF using |alpha_one| \newline */
  std::vector<double> cdf_one {};
  /* CDF using |alpha_two| \newline */
  std::vector<size_t> index {} ;
  std::vector<double> cdf_two {};
  std::vector<double> vEall {};
  std::vector<double> vE {}; 

    void freemem(){
      index.clear();
      std::vector<size_t>().swap( index);
    vE.clear();
    std::vector<double>().swap(vE);
    vEall.clear();
    std::vector<double>().swap(vEall);
    cdf_two.clear();
    std::vector<double>().swap(cdf_two);
    cdf_one.clear();
    std::vector<double>().swap(cdf_one); }

/* compute the mass function in Eq \eqref{PXN} \newline */
  double masspx( const double c_k, const double  c_alpha, const double c_cutoff )
  {
    return ( ( pow( 1.0/c_k, c_alpha) - pow( 1.0/(c_k+1.0), c_alpha) )/( 1.0 - pow( 1.0/( c_cutoff +1.0), c_alpha)) ); 
  }
};


@*1 {\bf the random number generator}. 
\label{SEC:rng}


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


@*1 {\bf compute the CDF}. 
\label{SEC:cdf}

Compute the cumulative density function for the distribution of juveniles Eq~\eqref{PXN} ; used for sampling number of juveniles using the inverse cdf method  in \S\~\ref{SEC:random_number_juveniles}.

@<cdf@>=@#
static void invcdf( M& self)
{
  self.cdf_one.clear();
  self.cdf_two.clear();

  self.cdf_one.push_back(0.);
  self.cdf_two.push_back(0.);

  double k =  1. ;
  while ( k <= self.psi){
   self.cdf_two.push_back(  self.cdf_two.back() +  self.masspx( k, self.alpha_two,  self.psi));
    self.cdf_one.push_back( self.cdf_one.back() +  self.masspx( k, self.alpha_one,  self.psi));
    ++k ; }

}


@*1 {\bf sample a random number of juveniles}.
\label{SEC:random_number_juveniles}

Sample a random number of juveniles using the inverse-cdf method,
i.e.\ drawing a random uniform and check where it lands, i.e.\ if $F$ denotes the cumulative density function  we compute
\begin{equation}
\label{eq:1}
j = \min\left\{ k \in \mathds{N} :  F(k) > u   \right\}
\end{equation}

The CDF is
computed in \S~\ref{SEC:cdf}. 

@<samplej@>=@#
static size_t samplexi( M& self, const  size_t one_two)
{

/* we sample at least one juvenile \newline */
  size_t j = 1;
  /* sample a random uniform \newline */
  const double u = gsl_rng_uniform( rngtype); ;
  if( one_two < 2 ){
  /* use CDF for |alpha_one| \newline */
  while( u > self.cdf_one[j] ){
    ++j;}}
    else {
      /* use CDF for |alpha_two| \newline */
    while( u > self.cdf_two[j] ){
    ++j;}
    }

  assert( j <= static_cast<size_t>( self.psi) ); 
  
  return (j) ;
}


@*1 {\bf comparison module for partial sorting}.
\label{SEC:comp}

Comparison module for partial sorting of juveniles given their  viability weight;  used in \S~\ref{SEC:nth}. 

@<comp@>=@#
static bool comp(const double a, const double b) 
{ 
    return (a < b); 
}


@*1 {\bf compute the $N$th  element }. 
\label{SEC:nth}

Compute the $N$th element using partial sorting into ascending order
using the comparison function \S~\ref{SEC:comp}. If $S_N > N$
juveniles, then returns the value $t$ so that
\begin{equation}
 \sum_{1 \le i \le S_N} \one{ t_i \le t  } = N
\end{equation}

@<nth@>=@#
static void nthelm( M& self )
{

  /* call from R as c(0, Rvec)  */
  std::nth_element( self.vEall.begin(), self.vEall.begin() + (self.N - 1), self.vEall.end(), comp); 

  self.Nth =  self.vEall[ self.N - 1] ;
}


@*1 {\bf sample pool of juveniles}. 
\label{SEC:samplepool}


Sample a random number of juveniles for a given subset of the current individuals using \S~\ref{SEC:random_number_juveniles}; returns the total number of juveniles. 

@<samplepooljuveniles@>=@#
static size_t samplepool( const size_t c_ninds,  const size_t c_one_two, M& self)
{
/* |c_ninds| is number of individuals of a given type; |c_one_two| is the indicator for using $\alpha_1$ or $\alpha_2$ \newline */
  assert( c_ninds > 0);
  assert( c_ninds <= self.N); 
  int j = static_cast<int>(c_ninds+1);
  size_t telja= 0;
  size_t SN = 0 ;
  
    while( --j > 0){
    /* counter for internal control \newline */
     ++telja ;
     /* |samplexi| in \S~\ref{SEC:random_number_juveniles} \newline */
     SN +=  samplexi(self, c_one_two); }
   assert( SN > 0 );
   assert( telja == c_ninds);

return( SN);
}



@*1 {\bf assign weight to juveniles}.
\label{SEC:assignweight}

Assign weight to all juveniles; the weight of the advantageous type is
one, of the wild type is $\exp(-s)$ where $s > 0$ is the strength of
selection. We record a random exponential with rate the corresponding weight.  

@<assignweight@>=@#
static void assignweight(M& self,  gsl_rng * r)
{
  /* assign weights to all juveniles; |SNW| is the number of juveniles of the advantageous type, |SN| the number of juveniles of the wild type  \newline */
  assert( self.SNW > 0);
  assert( self.SN > 0); 
   int i = static_cast<int>(self.SNW + 1);
   /* clear the containers with the weights \newline */
  self.vE.clear();
  /* |vEall| contains all the weights for  computing the $N$th element \newline */
  self.vEall.clear() ;
  self.vE.shrink_to_fit();
  std::vector<double>().swap( self.vE) ;
  self.vEall.shrink_to_fit();
  std::vector<double>().swap( self.vEall) ;
   while( --i > 0){
   /* the optimal genotype correspondingly trait value is  $g_0 = 0$, so weight is $W = 1$ \newline */
     self.vE.push_back(  gsl_ran_exponential( r,  1. ) );
     self.vEall.push_back( self.vE.back() ); }


   i = static_cast<int>( self.SN + 1);
    while( --i > 0){
    /* the wild type is denoted  1, so weight is $W = e^{-s}$ \newline */ 
       self.vEall.push_back( gsl_ran_exponential( r,  1./exp( - self.c_strength_selection ) ) );  }
}


@*1 {\bf count survivors of a bottleneck}. 
\label{SEC:survivingb}

Count number of individuals of the advantageous type surviving a
bottleneck $N_b$. We suppose the $y$ individuals of the advantageous type
are enumerated from $0$ to $y-1$, we shuffle the $N$ indexes and  count, with $p := \min( y, N_b )$,  
\begin{equation}
\label{eq:2}
y^\prime =  \sum_{i=1}^N \one{\sigma(i) < p}
\end{equation}
where $\sigma(i)$ is the shuffled index. 

@<surviving@>=@#
static void count_number_surviving_bottleneck(M& self)
{

  size_t newy = 0;
  const size_t reference_point = self.y < self.bottlenecksize ? self.y : self.bottlenecksize ; 
  std::random_shuffle( std::begin( self.index), std::end(self.index) ) ;
  for( const auto &i : self.index){
  /* check if an advantageous individual survives the bottleneck Eq~(\ref{eq:2}) \newline */
    newy += (i < reference_point ? 1 : 0); }

  /* set number of good type to count of surviving bottleneck \newline */
  self.y = (newy < self.bottlenecksize ? newy : self.N) ;
  self.Nprime = (newy <  self.bottlenecksize ? self.bottlenecksize - newy : 0)  ; 
}


@*1 {\bf surviving bottleneck by hypergeometric}. 
\label{sec:survivehypergeom}

The number surviving a bottleneck is a hypergeometric. 

@<Bottleneckypergeometric@>=@# 
static void surviving_bottleneckHypergeometric(M& self, gsl_rng *r)
{
 const unsigned  newy = gsl_ran_hypergeometric( r, self.y,  self.Nprime, self.bottlenecksize) ;
 
  /* set number of good type to count of surviving bottleneck \newline */
  self.y = (newy < self.bottlenecksize ? newy : self.N) ;
  /* |Nprime| is  number of wild type \newline */
  self.Nprime = (newy <  self.bottlenecksize ? self.bottlenecksize - newy : 0)  ; 
}


@*1 {\bf count surviving by weight}.
\label{sec:surviveweight}

Count the number of juveniles of the advantageous type  surviving a
selection by weight.

@<byweight@>=@#
static void count_number_surviving_assigning_weight(M& self)
{

  self.y = 0;
  /* |self.vE| is container with juveniles of the advantageous type
  \newline */
   for( const auto &t : self.vE){
    self.y  += (t <= self.Nth ? 1 : 0);}

 self.y =  (self.y <= self.N ? self.y : self.N ) ;
}


@*1 {\bf one step}.
\label{sec:onestep}

take one step through an excursion

@<onestep@>=@#
static void  onestep_after_bottleneck(M& self, gsl_rng * r)
{

  /* check if bottleneck occurs, and if it occurs sample surviving  \S~\ref{sec:survivehypergeom} \newline  */
  if( gsl_rng_uniform(r) < self.pbottle ){
    surviving_bottleneckHypergeometric(self, r) ;}
  if( self.y > 0){
    if( self.y < self.N ){
      /* need to sample juveniles \newline  */
       const size_t c_one_two = (gsl_rng_uniform(r) < self.epsilon_N ? 1 : 2) ;
       assert( self.y > 0);
       assert( self.Nprime > 0);
       /* sample juveniles with advantageous type \S~\ref{SEC:samplepool} \newline */
       self.SNW = samplepool( self.y, c_one_two, self);
       /* sample juveniles with wild type \newline */
       self.SN =   samplepool( self.Nprime, c_one_two, self);
	if( self.SNW + self.SN <= self.N ){
	  /* total number of juveniles does not exceed $N$ so all
          survive \newline */
	  self.y = self.SNW ;
	  self.Nprime = self.SN ; }
	else{
	  /* total number of juveniles exceeds $N$ so need to assign
          weights \S~\ref{SEC:assignweight} \newline */
	   assignweight( self, r);
           /* compute the $N$th smallest exponential \S~\ref{SEC:nth}
           \newline */ 
	   nthelm( self);
           /* \S~\ref{sec:surviveweight} \newline */
	   count_number_surviving_assigning_weight(self);
           assert( self.y <= self.N ) ;
	   self.Nprime = self.N - self.y ;
	}
    }
    else{
      assert( self.y == self.N); }
  }
  else{
    assert( self.y  < 1); }
}


@*1 {\bf trajectory }. 
\label{sec:trajectory}

@<traject@>=@#
static void trajectory ( M& self,  gsl_rng * r)
{
  std::vector<double> excursion_to_fixation {} ;
  double  timi = 0.0;
  int j = 0 ; ;
  std::vector<double> current_number_individuals {} ; 
  self.y = 1;
  self.Nprime = self.N - 1 ;
  while( ( self.y < self.N) && ( self.y > 0) ){
    current_number_individuals.push_back( self.Nprime + self.y); 
      excursion_to_fixation.push_back(self.y);
       timi +=  1.0;
       onestep_after_bottleneck(self, r); }

  printf("%d %g\n", self.y < self.N ? 0 : 1,  timi) ;

const std::string eskra = self.c_excursion_skra + std::to_string(self.numer) + ".txt" ; 
  if( self.y > 0){
   std::cout << eskra << '\n' ;
  std::ofstream f(eskra, std::ofstream::app);
  assert( f.is_open() ) ;
    for( auto const &y : excursion_to_fixation){
    f << y/current_number_individuals[j] << ' ';
    ++j;}
    f << "1\n";
    f.close();
  }
  excursion_to_fixation.clear();
  excursion_to_fixation.shrink_to_fit();
  std::vector<double>().swap( excursion_to_fixation);
}




@*1 {\bf the main module}. 
\label{SEC:main}


The |main| function

@C

@<includes@>@#
@<structM@>@#
@<gslrng@>@#
@<cdf@>@#
@<samplej@>@#
@<comp@>@#
@<nth@>@#
@<samplepooljuveniles@>@#
@<assignweight@>@#
@<surviving@>@#
@<Bottleneckypergeometric@>@#
@<byweight@>@#
@<onestep@>@#
@<traject@>@#

int main(int argc, char *argv[])
{
  M d {} ;
  
  setup_rng(  static_cast<unsigned long int>(atoi(argv[1])) );
  d.numer = atoi(argv[1]);

  invcdf( d);
  for( int i = 0 ; i < d.number_experiments ; ++i){
  printf("%d / %d : ", i, d.number_experiments ) ;
    trajectory( d, rngtype); }

/* clear the memory occupied by the containers \ref{SEC:structM} \newline */
d.freemem();
/* free the random number generator in \S~\ref{SEC:rng}  \newline */
  gsl_rng_free( rngtype ); 
  return 0 ;
}


@* {\bf example excursions}.
\label{sec:examples}


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen_}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=1.05$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.1$,
size  of bottleneck $10^{2}$, probability of a bottleneck in any given
generation  $0.01$}
\label{fig:examples}
\end{figure}

\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__B}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=1.05$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.1$}
\label{fig:exampleB}
\end{figure}


\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__C}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=1.05$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 0.1$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{3}$ experiments }
\label{fig:exampleC}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__D}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 0.1$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.01$; results from $10^{2}$ experiments }
\label{fig:exampleD}
\end{figure}


\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__E}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 0.1$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{2}$ experiments }
\label{fig:exampleE}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__F}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.5$,
size  of bottleneck $10^{2}$, probability of a bottleneck in any given
generation  $0.01$; results from $10^{2}$ experiments }
\label{fig:exampleF}
\end{figure}


\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__G}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{3}$ experiments  }
\label{fig:exampleG}
\end{figure}


\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__H}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{3}$ experiments  }
\label{fig:exampleH}
\end{figure}


\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__I}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.5$,
size  of bottleneck $10^{4}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{2}$ experiments}
\label{fig:exampleI}
\end{figure}

\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__J}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 1/N$,   selection strength $s=0.5$,
size  of bottleneck $10^{3}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{3}$ experiments }
\label{fig:exampleJ}
\end{figure}

\clearpage
\pagebreak
\newpage

\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__K}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 0.1$,   selection strength $s=0.5$,
size  of bottleneck $10^{3}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{3}$ experiments }
\label{fig:exampleK}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{excursion_skra__tablen__L}
\caption{Examples of excursions to fixation for $N=10^{6}$,
$\alpha_{1}=0.75$, $\alpha_{2} =3$, cutoff $\Psi_{N} = N$
Eq~\eqref{PXN}, $\varepsilon_{N} = 0.1$,   selection strength $s=0.5$,
size  of bottleneck $10^{2}$, probability of a bottleneck in any given
generation  $0.1$; results from $10^{3}$ experiments }
\label{fig:exampleL}
\end{figure}

\clearpage
\pagebreak
\newpage


@* {\bf conclusion}. 
\label{sec:concl}


We are interested in understanding how  random sweepstakes and
randomly occurring  bottlenecks  affect  viability selection in a
haploid population.  


In the absence of selection and bottlenecks the transition of $Y_{t}$
is given by a hypergeometric, for $0 \le x \le N$
\begin{equation}
\label{eq:5}
\prb{Y_{t+1} = x : Y_{t} = y} = \EE{ \prb{Y_{t+1} = x : Y_{t} = y,  S_{y}, S_{N - y}}} =\EE{ \frac{ \binom{S_{y}}{x} \binom{S_{N - y}}{N-x}  }{ \binom{S_{N}}{N}  }}
\end{equation}
where $S_{y}$ is the random number of juveniles produced by $y$
parents of the advantageous type, and $S_{N - y}$ is the random
number of juveniles produced by $N-y$ parents of the wild type and $S_{N} =
S_{Y_{t}} + S_{N - Y_{t}}$. Given our model of random sweepstakes
Eq~\eqref{eq:3} and Eq~\eqref{PXN}   we have  $S_{N} \ge N$ almost surely.  


In the presence of bottlenecks let $N_{t}$ denote the  population size
 at time $t$, i.e.\ with $B$ the bottleneck size $B \le  N_{t } \le N$
 and $1 < B < N$.  Then, with $S_{N_{t}} =
 S_{y} + S_{N_{t} - y}$
 \begin{equation}
 \label{eq:6}
 \prb{Y_{t+1} = x : Y_{t} = y} = \EE{\one{S_{y} = x, S_{N_{t}} \le N}} +  \EE{ \frac{ \binom{S_{y}}{x} \binom{S_{N_{t} - y}}{N-x}  }{ \binom{S_{N_{t}}}{N}  }\one{S_{N_{t}} > N} }
\end{equation}
Including selection corresponds to  weighted sampling where the
probability of sampling a juvenile with weight one is
\begin{equation}
\label{eq:7}
  \frac{S_{y} }{S_{y} +  e^{-s}S_{N_{t} - y} }
\end{equation}
when $y$ parents are of type with weight one.   

A mathematical investigation of our model will have to  be postponed
for future study.






@* {\bf references}.
\label{sec:refs}

\bibliographystyle{plain}
\bibliography{/home/bjarki/verk/master_bibfile/refs}



@
\end{document}