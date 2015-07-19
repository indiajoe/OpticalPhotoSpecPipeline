#Optical Imaging and Spectroscopy Pipeline
=========================================

This is a comprehensive data reduction pipeline for optical instruments.
Aim of this pipeline is to provide users a single pipeline interface for reducing multiple instrument data.

Currently, this pipeline supports following instruments
+ HFOSC @HCT, India : (Both Imaging and Spectroscopy data)
+ IFOSC @IGO, India : (Both Imaging and Spectroscopy data)


Main Documentation : https://github.com/indiajoe/OpticalPhotoSpecPipeline/wiki

The pipeline's user interface is very similar to [TIRSPEC](http://indiajoe.github.io/TIRSPEC/) (NIR instrument @HCT)
Hence, for some detailed instruction on how to run the pipeline, see its wiki page also : https://github.com/indiajoe/TIRSPEC/wiki

Feedbacks in the form of suggestions and pull requests are more than welcome.

## Why a single pipeline to reduce all data?
+ With commissioning of each new telescope and instrument, we are facing an exponential increase in astronomy data. Current strategy of reducing each frame manually is not scalable to future. Pipelines are the only way to tackle the data explosion.
+ With more facilities coming up, astronomers have access to data from different instruments. It is a painful job to learn different pipelines for different instruments to fundamentally do the same job. 
 
## Philosophy
+ This is *not* a quick look tool. The aim of this pipeline is to make the reduction as easy and streamlined as possible without compromising on the quality of reduction. We hope to do better than naive manual reduction by collectively putting together best practices of reduction.
+ Source code will always remain free and open. Any science done with closed source software/procedure is simply not science, full stop.

## Development Plan
+ As new instruments are added, we will have to progressively abstractify the code to accommodate them all. So the road-map is a slow an steady progress towards the most general framework.
