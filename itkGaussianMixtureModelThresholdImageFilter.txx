/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianMixtureModelThresholdImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2006/03/15 01:57:09 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkGaussianMixtureModelThresholdImageFilter_txx
#define _itkGaussianMixtureModelThresholdImageFilter_txx

#include "itkGaussianMixtureModelThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkGaussianMixtureModelThresholdCalculator.h"

namespace itk {

template<class TInputImage, class TOutputImage>
GaussianMixtureModelThresholdImageFilter<TInputImage, TOutputImage>
::GaussianMixtureModelThresholdImageFilter()
{
  m_OutsideValue   = NumericTraits<OutputPixelType>::Zero;
  m_InsideValue    = NumericTraits<OutputPixelType>::max();
  m_Threshold      = NumericTraits<InputPixelType>::Zero;
  m_NumberOfHistogramBins = 128;
}

template<class TInputImage, class TOutputImage>
void
GaussianMixtureModelThresholdImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typename ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  // Create a histogram of the image intensities
  typename HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
  histogramGenerator->SetInput(  this->GetInput()  );
  histogramGenerator->SetNumberOfBins( m_NumberOfHistogramBins );
  histogramGenerator->Compute();

  // Compute the multiple Otsu Thresholds for the input image
  typename CalculatorType::Pointer thresholdCalculator = CalculatorType::New();
  thresholdCalculator->SetInputHistogram( histogramGenerator->GetOutput() );
  thresholdCalculator->Update();

  m_Threshold = thresholdCalculator->GetOutput();

  typename BinaryThresholdImageFilter<TInputImage,TOutputImage>::Pointer threshold = 
    BinaryThresholdImageFilter<TInputImage,TOutputImage>::New();;
  
  progress->RegisterInternalFilter(threshold,.5f);
  threshold->GraftOutput (this->GetOutput());
  threshold->SetInput (this->GetInput());
  threshold->SetLowerThreshold( m_Threshold );
  threshold->SetInsideValue (m_InsideValue);
  threshold->SetOutsideValue (m_OutsideValue);
  threshold->Update();

  this->GraftOutput(threshold->GetOutput());
}

template<class TInputImage, class TOutputImage>
void
GaussianMixtureModelThresholdImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  const_cast<TInputImage *>(this->GetInput())->SetRequestedRegionToLargestPossibleRegion();
}

template<class TInputImage, class TOutputImage>
void 
GaussianMixtureModelThresholdImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Threshold: " << static_cast<typename NumericTraits<InputPixelType>::PrintType>(m_Threshold) << std::endl;
}


}// end namespace itk
#endif
