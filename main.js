  // Load MODIS image collections
  var modis = ee.ImageCollection('MODIS/061/MOD13Q1').filterBounds(geometry);
  var modistem = ee.ImageCollection("MODIS/061/MOD11A1").filterBounds(geometry);
  
  // Define date ranges for each year
  var startYear = 2012;
  var endYear = 2022;
  var dateRanges = ee.List.sequence(startYear, endYear).map(function(year) {
    return ee.Date.fromYMD(year, 1, 1).format('YYYY'); // Only format as year
  });
  
  // Cloud masking functions
  function filterClouds(image) {
    var qa = image.select('SummaryQA');
    var cloudState = qa.bitwiseAnd(3);
    var clearMask = cloudState.eq(0);
    return image.updateMask(clearMask);
  }
  
  function filterCloudsLST(image) {
    var qa = image.select('QC_Day');
    var cloudState = qa.bitwiseAnd(3);
    var clearMask = cloudState.eq(0);
    return image.updateMask(clearMask);
  }
  
  // Calculate 16-day EVI composite for a given year within geometry
  function get16DayEVIComposite(startDate) {
    var daysInterval = 16;
    var dates = ee.List.sequence(0, 365, daysInterval).map(function(dayOffset) {
      return ee.Date(startDate).advance(dayOffset, 'day');
    });

  var composite = ee.ImageCollection(dates.map(function(date) {
    var start = ee.Date(date);
    var end = start.advance(daysInterval, 'day');
    var evi = modis.filterDate(start, end)
      .map(filterClouds)
      .median()
      .select('EVI')
      .multiply(0.0001);
    return evi.set('system:time_start', start.millis());
  }));

    var yearlyComposite = composite.median().clip(geometry);
    return yearlyComposite;
  }
  
  // Calculate 16-day LST composite for a given year within geometry
  function get16DayLSTComposite(startDate) {
    var daysInterval = 16;
    var dates = ee.List.sequence(0, 365, daysInterval).map(function(dayOffset) {
      return ee.Date(startDate).advance(dayOffset, 'day');
    });
  
    var composite = ee.ImageCollection(dates.map(function(date) {
      var start = ee.Date(date);
      var end = start.advance(daysInterval, 'day');
      var lst = modistem.filterDate(start, end)
        .map(filterCloudsLST)
        .median()
        .select('LST_Day_1km')
        .multiply(0.02)
        .subtract(273.15);
      return lst.set('system:time_start', start.millis());
    }));
  
    var yearlyComposite = composite.median().clip(geometry);
    return yearlyComposite;
  }
  
  // Load MODIS Land Cover image collection (MCD12Q1) for masking forest areas
  var landCover = ee.ImageCollection('MODIS/061/MCD12Q1')
                  .filter(ee.Filter.date('2020-01-01', '2020-12-31'))
                  .first()
                  .select('LC_Type1');
  
  var forestTypes = ee.List([1, 2, 3, 4, 5]);
  var forestMask = landCover.remap({
    from: forestTypes,
    to: ee.List.repeat(1, forestTypes.length())
  }).eq(1);
  
  // Calculate Time Period 2 (entire span) within geometry
  var period2EVI = modis.filterDate('2012-01-01', '2022-12-31')
    .map(filterClouds)
    .median()
    .select(['EVI'])
    .multiply(0.0001)
    .clip(geometry);
  
  var period2LST = modistem.filterDate('2012-01-01', '2022-12-31')
    .map(filterCloudsLST)
    .median()
    .select(['LST_Day_1km'])
    .multiply(0.02)
    .subtract(273.15)
    .clip(geometry);
  
  // Initialize ImageCollections for EVI and LST differences
  var eviDifferenceCollection = ee.ImageCollection([]);
  var lstDifferenceCollection = ee.ImageCollection([]);
  
  // Map function to create collections of EVI and LST differences for each year
  dateRanges.evaluate(function(ranges) {
    ranges.forEach(function(startDate) {
      var startDateEE = ee.Date(startDate);
  
      // Calculate differences for EVI and LST for the current year
      var eviDifference = get16DayEVIComposite(startDateEE).subtract(period2EVI)
        .set('system:time_start', startDateEE.millis());
      var lstDifference = get16DayLSTComposite(startDateEE).subtract(period2LST)
        .set('system:time_start', startDateEE.millis());
  
      // Add the differences to their respective collections
      eviDifferenceCollection = eviDifferenceCollection.merge(ee.ImageCollection([eviDifference]));
      lstDifferenceCollection = lstDifferenceCollection.merge(ee.ImageCollection([lstDifference]));
    });
  
    // Calculate standard deviation for each difference collection
    var eviStdDev = eviDifferenceCollection.reduce(ee.Reducer.stdDev());
    var lstStdDev = lstDifferenceCollection.reduce(ee.Reducer.stdDev());
  
    // Calculate Z-Score for each difference collection
    var eviZScore = eviDifferenceCollection.map(function(image) {
      return image.divide(eviStdDev)
        .set('system:time_start', image.get('system:time_start'));
    });
  
    var lstZScore = lstDifferenceCollection.map(function(image) {
      return image.divide(lstStdDev)
        .set('system:time_start', image.get('system:time_start'));
  });
    
  // Rename bands in both Z-Score collections to have a common name
  var renamedEviZScore = eviZScore.map(function(image) {
    return image.rename('ZScore');
  });
  
  var renamedLstZScore = lstZScore.map(function(image) {
    return image.rename('ZScore');
  });
     
  // Calculate annual mean Z-Score for each collection
  function calculateAnnualMeanZScore(year, zScoreCollection, bandName) {
    var startDate = ee.Date.fromYMD(year, 1, 1);
    var endDate = ee.Date.fromYMD(year, 12, 31);
  
    // Calculate mean Z-Score for the given year
    var annualComposite = zScoreCollection.filterDate(startDate, endDate)
      .mean()
      .reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geometry,
        scale: 1000,
        maxPixels: 1e9
      });

  // Create a dictionary to hold the values
  var featureProps = {
    'year': ee.Number(year)
  };

  // Get the Z-score value, and ensure it's not null
  var zScoreValue = ee.Number(annualComposite.get('ZScore'));
  featureProps[bandName] = ee.Algorithms.If(zScoreValue, zScoreValue, null);

  // Return the feature with the constructed dictionary
  return ee.Feature(null, featureProps);
  }
  
  // Create a feature collection for each Z-Score time series
  var years = ee.List.sequence(startYear, endYear);
  var eviZScoreCollectionAnnual = ee.FeatureCollection(
    years.map(function(year) {
      return calculateAnnualMeanZScore(year, renamedEviZScore, 'EVI_ZScore');
    })
  );
  
  var lstZScoreCollectionAnnual = ee.FeatureCollection(
    years.map(function(year) {
      return calculateAnnualMeanZScore(year, renamedLstZScore, 'LST_ZScore');
    })
  );
  
  // Print intermediate results to diagnose
  print('EVI Z-Score Collection Annual:', eviZScoreCollectionAnnual);
  print('LST Z-Score Collection Annual:', lstZScoreCollectionAnnual);
  
  // Merge the collections into one for charting
  var mergedCollection = eviZScoreCollectionAnnual.map(function(feature) {
    var correspondingLST = lstZScoreCollectionAnnual.filter(ee.Filter.eq('year', feature.get('year'))).first();
    return feature.set('LST_ZScore', correspondingLST.get('LST_ZScore'));
  });
  
  // Plot the results
  var zScoreChart = ui.Chart.feature.byFeature(mergedCollection, 'year', ['EVI_ZScore', 'LST_ZScore'])
    .setChartType('LineChart')
    .setOptions({
      title: 'Annual Median Z-Scores for EVI and LST (2012-2022)',
      hAxis: {title: 'Year'},
      vAxis: {title: 'Mean Z-Score'},
      series: {
        0: {color: 'green', lineWidth: 2, pointSize: 5}, // EVI Z-Score
        1: {color: 'red', lineWidth: 2, pointSize: 5}    // LST Z-Score
      },
      legend: {position: 'bottom'}
  });
    
    // Display the chart
    print(zScoreChart);
    
  // Visualization parameters for Z-Score display
  var zScoreVisParams = {
    min: -2,
    max: 2,
    palette: ['blue', 'white', 'red']
  };
  
  // Initialize collections for Z-Score images
  var eviZScoreCollection = ee.ImageCollection([]);
  var lstZScoreCollection = ee.ImageCollection([]);
  
  // Calculate Z-Scores for each year and store in collections
  dateRanges.evaluate(function(ranges) {
    ranges.forEach(function(startDate) {
      var startDateEE = ee.Date(startDate);
  
      // Calculate differences for EVI and LST for the current year
      var eviDifference = get16DayEVIComposite(startDateEE).subtract(period2EVI)
        .set('system:time_start', startDateEE.millis());
      var lstDifference = get16DayLSTComposite(startDateEE).subtract(period2LST)
        .set('system:time_start', startDateEE.millis());
  
      // Calculate Z-Scores
      var eviZScoreImage = eviDifference.divide(eviStdDev).rename('EVI_ZScore');
      var lstZScoreImage = lstDifference.divide(lstStdDev).rename('LST_ZScore');
  
      // Add Z-Score images to collections
      eviZScoreCollection = eviZScoreCollection.merge(ee.ImageCollection([eviZScoreImage]));
      lstZScoreCollection = lstZScoreCollection.merge(ee.ImageCollection([lstZScoreImage]));
  

  });

  // Function to combine EVI and LST Z-Scores by iterating through each list
  var combinedZScores = ee.ImageCollection(ee.List.sequence(0, eviZScoreCollection.size().subtract(1)).map(function(index) {
    var eviZScore = ee.Image(eviZScoreCollection.toList(eviZScoreCollection.size()).get(index));
    var lstZScore = ee.Image(lstZScoreCollection.toList(lstZScoreCollection.size()).get(index));
    return eviZScore.multiply(lstZScore); // Multiply corresponding EVI and LST Z-Score images
  }));

  // EVI x LST Z-Scores, then sum
  var XX = combinedZScores.sum(); // Sum over the years

  // EVI Z-Scores squared, then sum
  var sevi2 = eviZScoreCollection.map(function(image) {
    return image.pow(2); // Square EVI Z-Score
  }).sum(); // Sum over the years

  // LST Z-Scores squared, then sum
  var slst2 = lstZScoreCollection.map(function(image) {
    return image.pow(2); // Square LST Z-Score
  }).sum(); // Sum over the years

  // Calculate the square root of (sevi2 * slst2)
  var YY = sevi2.multiply(slst2).sqrt();

  // Final calculation: XX / YY
  var final = XX.divide(YY);

  // Apply forest mask to the final result
  var maskedFinal = final.updateMask(forestMask);

  // Print or visualize the result
  print(maskedFinal);

  // Add the result to the map for visualization
  Map.addLayer(maskedFinal, {min: -1, max: 1, palette: ['blue', 'white', 'red']}, 'Correlation (Forest Masked)');

 
  });
  
  
  // Display a time-series chart for EVI Difference with time
  var eviChart = ui.Chart.image.series({
    imageCollection: eviDifferenceCollection,
    region: geometry, // Define your region of interest
    reducer: ee.Reducer.mean(), // Use mean value for pixel size
    scale: 1000, // Pixel size/value
    xProperty: 'system:time_start' // This property specifies the x-axis as the time property
  }).setOptions({
    series: {
      0: {color: 'blue'} // Set line color for EVI Difference
    },
    lineWidth: 1, // Set this to a non-zero value to connect points with lines
    pointSize: 1, // Adjust the point size
    title: 'EVI Difference Time Series',
    interpolateNulls: false,
    vAxis: {title: 'EVI Difference'},
    hAxis: {title: 'Date', format: 'YYYY'}
  });

  print(eviChart); // Display the chart in the console

  // Display a time-series chart for LST Difference with time
  var lstChart = ui.Chart.image.series({
    imageCollection: lstDifferenceCollection,
    region: geometry, // Define your region of interest
    reducer: ee.Reducer.mean(), // Use mean value for pixel size
    scale: 1000, // Pixel size/value
    xProperty: 'system:time_start' // This property specifies the x-axis as the time property
  }).setOptions({
    series: {
      0: {color: 'red'} // Set line color for LST Difference
    },
    lineWidth: 1, // Set this to a non-zero value to connect points with lines
    pointSize: 1, // Adjust the point size
    title: 'LST Difference Time Series',
    interpolateNulls: false,
    vAxis: {title: 'LST Difference (°C)'},
    hAxis: {title: 'Date', format: 'YYYY'}
  });

  print(lstChart); // Display the chart in the console
 
  
 
  // // Export the masked final result as an image 
  // Export.image.toDrive({
  //   image: maskedFinal,
  //   description: 'Pearson_Correlation',
  //   folder: 'EarthEngineExports',
  //   fileNamePrefix: 'Pearson_Correlation',
  //   region: geometry,
  //   scale: 500,
  //   maxPixels: 1e13,
  //   crs: 'EPSG:4326'
  // });

  // Export.image.toDrive({
  //   image: maskedFinal,
  //   description: 'Correlation_Forest_Masked',
  //   folder: 'EarthEngineExports',
  //   fileNamePrefix: 'Correlation_Forest_Masked',
  //   region: geometry,
  //   scale: 500,
  //   maxPixels: 1e13,
  //   crs: 'EPSG:4326'
  //   });
   
  //   // Export EVI difference for the current year 
  //   Export.image.toDrive({
  //     image: eviDifference,
  //     description: 'EVI_Difference_' + startDate, // Name file with year
  //     folder: 'EarthEngineExports',
  //     fileNamePrefix: 'EVI_Difference_' + startDate,
  //     region: geometry,
  //     scale: 500,
  //     maxPixels: 1e13,
  //     crs: 'EPSG:4326'
  //   });
   
  //   // Export LST difference for the current year
  //   Export.image.toDrive({
  //     image: lstDifference,
  //     description: 'LST_Difference_' + startDate, // Name file with year
  //     folder: 'EarthEngineExports',
  //     fileNamePrefix: 'LST_Difference_' + startDate,
  //     region: geometry,
  //     scale: 500,
  //     maxPixels: 1e13,
  //     crs: 'EPSG:4326'
  //   });
   
        // Export EVI Z-Score Image
  //    Export.image.toDrive({
  //      image: eviZScoreImage.clip(geometry),
  //      description: 'EVI_ZScore_' + startDate,
  //      folder: 'MODIS_ZScore_Exports', // Folder name in Google Drive
  //      fileNamePrefix: 'EVI_ZScore_' + startDate,
  //      region: geometry,
  //      scale: 1000,
  //      maxPixels: 1e9
  //    });
  
      // Export LST Z-Score Image
  //    Export.image.toDrive({
  //      image: lstZScoreImage.clip(geometry),
  //      description: 'LST_ZScore_' + startDate,
  //      folder: 'MODIS_ZScore_Exports', // Folder name in Google Drive
  //      fileNamePrefix: 'LST_ZScore_' + startDate,
  //      region: geometry,
  //      scale: 1000,
  //      maxPixels: 1e9
  //    });
   

});

Map.centerObject(geometry, 6);

 