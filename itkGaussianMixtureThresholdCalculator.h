/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianMixtureThresholdCalculator.h,v $
  Language:  C++
  Date:      $Date: 2005/11/25 15:50:35 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianMixtureThresholdCalculator_h
#define __itkGaussianMixtureThresholdCalculator_h

#include "itkHistogramAlgorithmBase.h"
#include "itkHistogram.h"

namespace itk
{

/** \class GaussianMixtureThresholdCalculator
 * \brief Computes Otsu's thresholds for a histogram.
 * 
 * You plug in the target histogram using SetInputHistogram method and 
 * specify the number of thresholds you want to be computed. Then call
 * the GenerateData method to run the alogithm. 
 *
 * The thresholds are computed so that the between-class variance is
 * maximized.
 *
 * \ingroup Calculators
 */

template< class TInputHistogram >
class GaussianMixtureThresholdCalculator :
      public HistogramAlgorithmBase< TInputHistogram >
{
public:
  /**Standard class typedefs. */
  typedef GaussianMixtureThresholdCalculator Self;
  typedef HistogramAlgorithmBase<TInputHistogram> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef typename TInputHistogram::MeasurementType MeasurementType;
  typedef typename TInputHistogram::FrequencyType FrequencyType;

  typedef typename NumericTraits<MeasurementType>::RealType MeanType;
  typedef typename NumericTraits<MeasurementType>::RealType VarianceType;

  typedef std::vector<MeanType> MeanVectorType;
  typedef std::vector<FrequencyType> FrequencyVectorType;

  typedef typename TInputHistogram::InstanceIdentifier InstanceIdentifierType;
  typedef std::vector<InstanceIdentifierType> InstanceIdentifierVectorType;

  /**Standard Macros */
  itkTypeMacro(GaussianMixtureThresholdCalculator, HistogramAlgorithmsBase);
  itkNewMacro(Self) ;

  /** Returns the threshold */
  itkGetConstMacro( Output, MeasurementType );

protected:
  GaussianMixtureThresholdCalculator() {};
  virtual ~GaussianMixtureThresholdCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Calculates the thresholds and save them */
  void GenerateData() ;

private:
  /** Internal thresholds storage */
  MeasurementType m_Output ;

} ; // end of class

} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianMixtureThresholdCalculator.txx"
#endif

#endif
