/****************************************************************************
*   vtkRegistrationQuantification.cxx
*
*   Created by:     Michael Kuczynski
*   Created on:     16/04/2019
*   Version:        1.0
*   Description:    
****************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include <vtkMath.h>
#include <vtkSmartPointer.h>

#include <vtkImageData.h>
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>

#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>

#include <vtkExtractVOI.h>
#include <vtkImageCast.h>

#include <vtkImageGaussianSmooth.h>
#include <vtkMarchingCubes.h>
#include <vtkImageThreshold.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkDecimatePro.h>
#include <vtkImageReslice.h>
#include <vtkImageChangeInformation.h>

#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>

#include <vtkImageMapper.h>
#include <vtkImageMapper3D.h>
#include <vtkActor.h>
#include <vtkActor2D.h>

#include <vtkProperty.h>
#include <vtkProperty2D.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>

static double ssd = 0;

static double sourceEntropy = 0;
static double targetEntropy = 0;
static double jointEntropy = 0;

void normalizeImage( vtkSmartPointer<vtkImageData> sourceImage );
void mutualInformation( vtkSmartPointer<vtkImageData> sourceImage, vtkSmartPointer<vtkImageData> targetImage, vtkSmartPointer<vtkImageData> outputImage );
void sumSquaredDifference( vtkSmartPointer<vtkImageData> sourceImage, vtkSmartPointer<vtkImageData> targetImage, vtkSmartPointer<vtkImageData> outputImage );
void correlationCoefficient( vtkSmartPointer<vtkImageData> sourceImage, vtkSmartPointer<vtkImageData> targetImage, vtkSmartPointer<vtkImageData> outputImage );

int main(int argc, char* argv[])
{
    // Path to MHA files
    std::string input_1 = argv[1];  // source image that has been transformed to the target image domain
    std::string input_2 = argv[2];  // target image

    vtkSmartPointer<vtkMetaImageReader> niftiReader = vtkSmartPointer<vtkMetaImageReader>::New();

    vtkSmartPointer<vtkImageData> sourceImage = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkImageData> targetImage = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkImageData> ssdImage    = vtkSmartPointer<vtkImageData>::New();

    niftiReader->SetFileName( input_1.c_str() );
    niftiReader->Update();
    sourceImage->DeepCopy(niftiReader->GetOutput());

    niftiReader->SetFileName( input_2.c_str() );
    niftiReader->Update();
    targetImage->DeepCopy(niftiReader->GetOutput());

    mutualInformation( sourceImage, targetImage, ssdImage );

//     // Create the LUT for overlaying the transformed source image on the target image
//     vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
//     lookupTable->SetNumberOfTableValues( 10 );
//     lookupTable->SetTableRange( 0.0, 11.0 );
//     // lookupTable->SetTableValue( 0.0, 0.0, 0.0, 0.0, 0.0 );
//     // lookupTable->SetTableValue( 1.0, 1.0, 0.0, 0.0, 1.0 );
//     lookupTable->Build();

//     vtkSmartPointer<vtkImageMapToColors> mapTransparency = vtkSmartPointer<vtkImageMapToColors>::New();
//     mapTransparency->SetLookupTable( lookupTable );
//     mapTransparency->PassAlphaToOutputOn();
//     mapTransparency->SetInputData( sourceImage );

//     // Setup mappers
//     vtkSmartPointer<vtkImageMapper> sourceMapper = vtkSmartPointer<vtkImageMapper>::New();
//     sourceMapper->SetInputConnection( mapTransparency->GetOutputPort() );
//     sourceMapper->SetZSlice( 150 );
//     sourceMapper->SetColorWindow( 1 );
//     sourceMapper->SetColorLevel( 1 );

//     vtkSmartPointer<vtkImageMapper> targetMapper = vtkSmartPointer<vtkImageMapper>::New();
//     targetMapper->SetInputData( targetImage );
//     targetMapper->SetZSlice( 150 );
//     targetMapper->SetColorWindow( 8650 );
//     targetMapper->SetColorLevel( 3295 );

//     vtkSmartPointer<vtkImageMapper> ssdMapper = vtkSmartPointer<vtkImageMapper>::New();
//     ssdMapper->SetInputData( ssdImage );
//     ssdMapper->SetZSlice( 150 );
//     ssdMapper->SetColorWindow( 100000 );
//     ssdMapper->SetColorLevel( 60000 );

//     // Setup actors
//     vtkSmartPointer<vtkActor2D> sourceActor = vtkSmartPointer<vtkActor2D>::New();
//     sourceActor->SetMapper( sourceMapper );
//     sourceActor->GetProperty()->SetOpacity( 0.5 );

//     vtkSmartPointer<vtkActor2D> targetActor = vtkSmartPointer<vtkActor2D>::New();
//     targetActor->SetMapper( targetMapper );
//     targetActor->GetProperty()->SetOpacity( 0.5 );

//     vtkSmartPointer<vtkActor2D> ssdActor = vtkSmartPointer<vtkActor2D>::New();
//     ssdActor->SetMapper( ssdMapper );

//     // Setup the renderer and render window
//     vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//     // renderer->AddActor( sourceActor );
//     // renderer->AddActor( targetActor );
//     renderer->AddActor( ssdActor );
//     renderer->ResetCamera();

//     vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//     renderWindow->AddRenderer( renderer );

//     vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//     vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();

//     renderWindowInteractor->SetInteractorStyle( style );

//     renderWindowInteractor->SetRenderWindow( renderWindow );
//     renderWindowInteractor->Initialize();

//     renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}

void sumSquaredDifference( vtkSmartPointer<vtkImageData> sourceImage, vtkSmartPointer<vtkImageData> targetImage, vtkSmartPointer<vtkImageData> outputImage )
{
    int sourceDimX = sourceImage->GetDimensions()[0];
    int sourceDimY = sourceImage->GetDimensions()[1];
    int sourceDimZ = sourceImage->GetDimensions()[2];

    outputImage->SetDimensions( sourceDimX, sourceDimY, sourceDimZ );
    outputImage->AllocateScalars(VTK_FLOAT, 1);

    int voxelCount = 0;

    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                float sourceVoxel = sourceImage->GetScalarComponentAsFloat(x, y, z, 0);
                float targetVoxel = targetImage->GetScalarComponentAsFloat(x, y, z, 0);
                float diff = std::abs( sourceVoxel - targetVoxel );
                float sqr = diff * diff;

                outputImage->SetScalarComponentFromFloat(x, y, z, 0, sqr);

                ssd += sqr;
                voxelCount++;
            }
        }
    }

    ssd /= voxelCount;

    std::cout << "SSD is: " << ssd << std::endl;

    // Write out the SSD image
    std::string filePath = "D:\\Git\\VTK_registration_quantification\\img\\ssdImage.nii";

    vtkSmartPointer<vtkNIFTIImageWriter > writer = vtkSmartPointer<vtkNIFTIImageWriter>::New();
    writer->SetInputData( outputImage );
    writer->SetFileName( filePath.c_str() );
    writer->Write();
}

void mutualInformation( vtkSmartPointer<vtkImageData> sourceImage, vtkSmartPointer<vtkImageData> targetImage, vtkSmartPointer<vtkImageData> outputImage )
{
    int sourceDimX = sourceImage->GetDimensions()[0];
    int sourceDimY = sourceImage->GetDimensions()[1];
    int sourceDimZ = sourceImage->GetDimensions()[2];

    outputImage->SetDimensions( 256, 256, sourceDimZ );
    outputImage->AllocateScalars(VTK_FLOAT, 1);

    int sourceMin = -1565;
    int sourceMax = 10343;

    int targetMin = -1325;
    int targetMax = 9993;

    int newMin = 0;
    int newMax = 255;

    // Normalize the images first to be 8-bit (intensity ranges from [0, 255])
    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                // Source image normalization
                float sourceVoxel = sourceImage->GetScalarComponentAsFloat(x, y, z, 0);
                float scale = (sourceVoxel - sourceMin) / ( sourceMax - sourceMin );
                float newSourceVoxel = ( scale * ( newMax - newMin ) ) + newMin;
            
                sourceImage->SetScalarComponentFromFloat(x, y, z, 0, newSourceVoxel);

                // Target image normalization
                // float targetVoxel = targetImage->GetScalarComponentAsFloat(x, y, z, 0);

                // sourceImage->SetScalarComponentFromFloat(x, y, z, 0, sqr);
            }
        }
    }

    std::cout << "Calculating joint entropy...";

    // Now calculate the joint entropy
    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                int sourceVoxel = int( sourceImage->GetScalarComponentAsFloat(x, y, z, 0) );
                int targetVoxel = int( sourceImage->GetScalarComponentAsFloat(x, y, z, 0) );
                float currentValue = outputImage->GetScalarComponentAsFloat(sourceVoxel, targetVoxel, z, 0);
                outputImage->SetScalarComponentFromFloat(sourceVoxel, targetVoxel, z, 0, currentValue + 1);
            }
        }
    }

    // Write out the MI image
    std::string filePath = "D:\\Git\\VTK_registration_quantification\\img\\miImage.nii";

    vtkSmartPointer<vtkNIFTIImageWriter > writer = vtkSmartPointer<vtkNIFTIImageWriter>::New();
    writer->SetInputData( sourceImage );
    writer->SetFileName( filePath.c_str() );
    writer->Write();

}

void normalizeImage( vtkSmartPointer<vtkImageData> sourceImage )
{
    int sourceDimX = sourceImage->GetDimensions()[0];
    int sourceDimY = sourceImage->GetDimensions()[1];
    int sourceDimZ = sourceImage->GetDimensions()[2];

    int sourceMin = -1565;
    int sourceMax = 10343;

    int targetMin = -1325;
    int targetMax = 9993;

    int newMin = 0;
    int newMax = 255;

    // Normalize the images first to be 8-bit (intensity ranges from [0, 255])
    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                // Source image normalization
                float sourceVoxel = sourceImage->GetScalarComponentAsFloat(x, y, z, 0);
                float scale = (sourceVoxel - sourceMin) / ( sourceMax - sourceMin );
                float newSourceVoxel = ( scale * ( newMax - newMin ) ) + newMin;
            
                sourceImage->SetScalarComponentFromFloat(x, y, z, 0, newSourceVoxel);

                // Target image normalization
                // float targetVoxel = targetImage->GetScalarComponentAsFloat(x, y, z, 0);

                // sourceImage->SetScalarComponentFromFloat(x, y, z, 0, sqr);
            }
        }
    }
}

void correlationCoefficient( vtkSmartPointer<vtkImageData> sourceImage, vtkSmartPointer<vtkImageData> targetImage, vtkSmartPointer<vtkImageData> outputImage )
{
    int sourceDimX = sourceImage->GetDimensions()[0];
    int sourceDimY = sourceImage->GetDimensions()[1];
    int sourceDimZ = sourceImage->GetDimensions()[2];

    // Calculate the SSD between the images and write the output image for visualization
    float meanA = 0;
    float meanB = 0;
    float varA = 0;
    float varB = 0;
    float stdA = 0;
    float stdB = 0;
    int countA = 0;
    int countB = 0;

    // Calculate the signal mean and standard deviation of each image
    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                float tempMeanA = sourceImage->GetScalarComponentAsFloat(x, y, z, 0);
                float tempMeanB = targetImage->GetScalarComponentAsFloat(x, y, z, 0);

                if ( tempMeanA > 0 )
                {
                    meanA += tempMeanA;
                    countA++;
                }

                if ( tempMeanB > 0 )
                {
                    meanB += tempMeanB;
                    countB++;
                }
            }
        }
    }

    meanA /= countA;
    meanB /= countB;

    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                float tempA = sourceImage->GetScalarComponentAsFloat(x, y, z, 0);
                float tempB = targetImage->GetScalarComponentAsFloat(x, y, z, 0);

                if ( tempA > 0 )
                {
                    varA += std::pow( ( tempA - meanA ), 2 );
                }

                if ( tempB > 0 )
                {
                    varB += std::pow( ( tempB - meanB ), 2 );
                }
            }
        }
    }

    varA /= countA;
    varB /= countB;

    stdA = sqrt( varA );
    stdB = sqrt( varB );

    float cc = 0;

    for ( int z = 0; z < sourceDimZ; z++ )
    {
        for ( int y = 0; y < sourceDimY; y++ )
        {
            for ( int x = 0; x < sourceDimX; x++ )
            {
                float sourceVoxel = sourceImage->GetScalarComponentAsFloat(x, y, z, 0);
                float targetVoxel = targetImage->GetScalarComponentAsFloat(x, y, z, 0);

                float temp = ( sourceVoxel - meanA ) * ( targetVoxel - meanB );
                temp /= sqrt( std::pow( ( sourceVoxel - meanA ), 2 ) * std::pow( ( targetVoxel - meanB ), 2 ) );

                cc += temp;

                outputImage->SetScalarComponentFromFloat(x, y, z, 0, temp);
            }
        }
    }

    // std::cout << "Correlation coefficient = " << cc << std::endl;
    // std::cout << "Mean of source image = " << meanA << std::endl;
    // std::cout << "Mean of target image = " << meanB << std::endl;
    // std::cout << "Standard deviation of source image = " << stdA << std::endl;
    // std::cout << "Standard deviation of target image = " << stdB << std::endl;
}