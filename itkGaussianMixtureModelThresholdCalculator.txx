/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianMixtureModelThresholdCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2005/01/13 15:36:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkGaussianMixtureModelThresholdCalculator_txx
#define _itkGaussianMixtureModelThresholdCalculator_txx

#include "itkGaussianMixtureModelThresholdCalculator.h"

namespace itk
{


/*
 * Compute Otsu's thresholds
 */                    
template<class TInputHistogram>
void
GaussianMixtureModelThresholdCalculator<TInputHistogram>
::GenerateData()
{
  typename TInputHistogram::ConstPointer histogram = this->GetInputHistogram();

/*std::cout << std::endl<< histogram->GetTotalFrequency()<< std::endl<< std::endl;
std::cout.flush();*/
// std::cout << std::endl<< histogram->GetSize()<< std::endl<< std::endl;
// std::cout << std::endl<< histogram->GetSize(0)<< std::endl<< std::endl;


// TODO: as an improvement, the class could accept multi-dimensional histograms
  // and the user could specify the dimension to apply the algorithm to.
  if (histogram->GetSize().GetSizeDimension() != 1)
    {
    itkExceptionMacro(<<"Histogram must be 1-dimensional.");
    }


  // iterate over all the bins in the histogram, excepted first one and last one,
  // estimate the parameters of the gaussian separated by that bin, and find
  // the bin which give the minimum error. This bin will be used for the thresholds

  double minError = NumericTraits<double>::max();
  unsigned long thIdx = 0;
  double thmu1 = 0;
  double thmu2 = 0;
  double thsigma1 = 0;
  double thsigma2 = 0;
  double thmult1 = 0;
  double thmult2 = 0;

  for( unsigned long idx=1; idx<histogram->Size()-1; idx++)
    {
    double mu1 = 0;
    double mu2 = 0;
    double card1 = 0;
    double card2 = 0;
    double sigma1 = 0;
    double sigma2 = 0;
    double mult1 = 0;
    double mult2 = 0;


    // compute mu and card for the 2 classes

    for( unsigned long i=0; i<=idx; i++)
      {
      card1 += histogram->GetFrequency( i );
      mu1 += histogram->GetMeasurementVector( i )[0] * histogram->GetFrequency( i );

// std::cout << "freq: " << histogram->GetFrequency( i ) << " bin v: " << histogram->GetMeasurementVector( i )[0] << " mu1: " << mu1 << " card1: " << card1 << std::endl;


      }

    for( unsigned long i=idx+1; i<histogram->Size(); i++)
      {
      card2 += histogram->GetFrequency( i );
      mu2 += histogram->GetMeasurementVector( i )[0] * histogram->GetFrequency( i );
      }

    if( card1 == 0 || card2 ==0 )
      {
      // a class is empty, no need to do the rest
      break;
      }

    // complete mu computation
    mu1 /= card1;
    mu2 /= card2;

    // compute sigma

    for( unsigned long i=0; i<=idx; i++)
      {
      sigma1 += histogram->GetFrequency( i ) * vnl_math_sqr( histogram->GetMeasurementVector( i )[0] - mu1 );
      }

    for( unsigned long i=idx+1; i<histogram->Size(); i++)
      {
      sigma2 += histogram->GetFrequency( i ) * vnl_math_sqr( histogram->GetMeasurementVector( i )[0] - mu2 );
      }

    // complete sigma computation
    sigma1 /= card1;
    sigma2 /= card1;

    // mult
    typename TInputHistogram::MeasurementVectorType mv;
    mv[0] = mu1;
    unsigned long hi = histogram->GetIndex( mv )[0];
    mult1 = NumericTraits< double >::NonpositiveMin();
    for( int i=hi-1; i<=hi+1; i++)
      { mult1 = vnl_math_max( (double)histogram->GetFrequency( i ), mult1 ); }

    mv[0] = mu2;
    hi = histogram->GetIndex( mv )[0];
    mult2 = NumericTraits< double >::NonpositiveMin();
    for( int i=hi-1; i<=hi+1; i++)
      { mult2 = vnl_math_max( (double)histogram->GetFrequency( i ), mult2 ); }
    

    // ok, we have every thing needed to describe our to gaussians
    // now estimate if they fits well with the data
    double error = 0;
    for( unsigned long i=0; i<=idx; i++)
      {
      double gamma1 = mult1 * vcl_exp( -vnl_math_sqr( histogram->GetMeasurementVector( i )[0] - mu1 ) / ( 2*sigma1 ) );
      double gamma2 = mult2 * vcl_exp( -vnl_math_sqr( histogram->GetMeasurementVector( i )[0] - mu2 ) / ( 2*sigma2 ) );
      double gamma = gamma1 + gamma2;
      
      error += vnl_math_sqr( gamma - histogram->GetFrequency( i ) );
      }

    if( error < minError )
      {
      minError = error;
      thIdx = idx;
      thmu1 = mu1;
      thmu2 = mu2;
      thsigma1 = sigma1;
      thsigma2 = sigma2;
      thmult1 = mult1;
      thmult2 = mult2;
      }
// std::cout << "idx: " << idx << " mu1: " << mu1 << " mu2: " << mu2 << " error: " << error << std::endl;

    }


  typename TInputHistogram::MeasurementVectorType mv;
  mv[0] = thmu1;
  unsigned long i1 = histogram->GetIndex( mv )[0];
  mv[0] = thmu2;
  unsigned long i2 = histogram->GetIndex( mv )[0];

  unsigned long minIdx = i1;
  double min = NumericTraits< double >::max();
  for( unsigned long i=i1; i<=i2; i++ )
    {
    double gamma1 = thmult1 * vcl_exp( -vnl_math_sqr( histogram->GetMeasurementVector( i )[0] - thmu1 ) / ( 2*thsigma1 ) );
    double gamma2 = thmult2 * vcl_exp( -vnl_math_sqr( histogram->GetMeasurementVector( i )[0] - thmu2 ) / ( 2*thsigma2 ) );
    double diffGamma2 = vnl_math_sqr( gamma1 - gamma2 );
    if( diffGamma2 < min )
      {
      min = diffGamma2;
      minIdx = i;
      }
    }

   m_Output = histogram->GetMeasurementVector( minIdx )[0];

// std::cout << "thidx: " << thIdx << std::endl;
// std::cout << "direct th: " << histogram->GetMeasurementVector( thIdx )[0] << std::endl;
// std::cout << "min th: " << histogram->GetMeasurementVector( minIdx )[0] << std::endl;
// std::cout << "mu1: " << thmu1 << std::endl;
// std::cout << "mu2: " << thmu2 << std::endl;
// std::cout << "sig1: " << thsigma1 << std::endl;
// std::cout << "sig2: " << thsigma2 << std::endl;

std::cout.flush();

}

template<class TInputHistogram>
void
GaussianMixtureModelThresholdCalculator<TInputHistogram>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Output: " << static_cast<typename NumericTraits<MeasurementType>::PrintType>(m_Output) << std::endl;

}

} // end namespace itk

#endif
