'use strict';

// Get the data from the static information compiled into the HTML page from the server;
// Also, prevent these properties to be changed 
// Object.freeze(appData); 

// Mapping of accepted energy units.
// Note: `from_meV_factor` is the factor to multiply to a meV value to obtain the required unit
var energyUnitsMapping = {
    meV: {
        label: "meV",
        from_meV_factor: 1.,
    },
    THz: {
        label: "THz",
        from_meV_factor: 0.241798926,
    },
    cm_minus_one: {
        label: "cm⁻¹",
        from_meV_factor: 8.06554399711,
    }
}

var app = new Vue({
    el: '#app',
    data: {
        maxLayers: 10,
        warningMessage: '',
        showSpinner: false,
        currentRequest: null,
        appData: appData,
        latexCode: "",
        modesFilterX: true,
        modesFilterY: true,
        modesFilterZ: true,
        forceConstantVariables: [],
        seriesWithData: {},
        forceConstantUnits: "1e19_N_over_m3",
        frequencyUnits: "cm_minus_one",

        seriesMetadata: Object.freeze([ // Initial series definition
            {
                internalName: 'IR_Y+Raman_Y+BackScattering_Y',
                name: 'Infrared + Raman active (visible in backscattering)',
                marker: {
                    symbol: 'mycrosscircle',
                    lineWidth: 2,
                    fillColor: 'rgba(255, 255, 255, 0)',
                    lineColor: 'rgb(223, 83, 83)'
                },
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_Y+Raman_Y+BackScattering_N',
                name: 'Infrared + Raman active (not visible in backscattering)',
                marker: {
                    symbol: 'mycrosscircle',
                    lineWidth: 2,
                    fillColor: 'rgba(255, 255, 255, 0)',
                    lineColor: 'rgb(55, 126, 184)'
                },
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_Y+Raman_N+BackScattering_Y',
                name: 'SHOULD_NOT_HAPPEN_Infrared active only (visible in backscattering)',
                marker: {
                    symbol: 'cross',
                    lineWidth: 2,
                    lineColor: 'rgb(223, 83, 83)'
                },
                color: 'rgb(223, 83, 83)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_Y+Raman_N+BackScattering_N',
                name: 'Infrared active only',
                marker: {
                    symbol: 'cross',
                    lineWidth: 2,
                    lineColor: 'rgb(55, 126, 184)'
                },
                color: 'rgb(55, 126, 184)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_Y+BackScattering_Y',
                name: 'Raman active only (visible in backscattering)',
                marker: {
                    symbol: 'mycircle',
                    lineWidth: 2,
                    fillColor: 'rgba(255, 255, 255, 0)',
                    lineColor: 'rgb(223, 83, 83)'
                },
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_Y+BackScattering_N',
                name: 'Raman active only (not visible in backscattering)',
                marker: {
                    symbol: 'mycircle',
                    lineWidth: 2,
                    fillColor: 'rgba(255, 255, 255, 0)',
                    lineColor: 'rgb(55, 126, 184)'
                },
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_N+BackScattering_Y',
                name: 'SHOULD_NOT_HAPPEN_Inactive (visible in backscattering)',
                marker: {
                    symbol: 'square'
                },
                color: 'rgb(223, 83, 83)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_N+BackScattering_N',
                name: 'Inactive',
                marker: {
                    symbol: 'square'
                },
                color: 'rgb(80, 80, 80)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            }
        ])
    },
    methods: {
        getLatexHTML: function(latexcode) {
            if (_.has(window, 'katex')) {
                return window.katex.renderToString(latexcode, {
                    throwOnError: false,
                    displayMode: true
                });
            }
            return '';
        },
        getUpdateDebouncedFunction: function() {
            this.showSpinner = true;
            return _.debounce(this.fetchAndUpdateData, 200); // Run the request only after 200ms, and don't rerun if new requests arrive (only run the last)
        },
        updateGraphDebounced: function() {
            // will be replaced later in 'mount', but needs to be there already
            // see https://vuejs.org/v2/guide/reactivity.html#Declaring-Reactive-Properties
        },
        fetchAndUpdateData: function() {
            try {
                this.currentRequest.abort();
            } catch(err) {
                if (err.name !== 'TypeError') {
                    throw err;
                }
            }
            this.showSpinner = true;
            this.warningMessage = '';
            $.ajaxSetup({
                contentType: "application/json"
              });
            var vueApp = this;

            var dataToSend = {
                forceConstantParams: _.object(
                    _.map(vueApp.appData.forceConstants.variables, function(variable) {return variable['name'];}),
                    _.map(vueApp.appData.forceConstants.variables, function(variable) {return variable['value'];})
                ),
                matrices: vueApp.appData.forceConstants.matrices,
                forceConstantUnits: vueApp.forceConstantUnits,
                pointgroupEven: vueApp.appData.pointgroupEven,
                pointgroupOdd: vueApp.appData.pointgroupOdd,
                uniqueAxisTransformationEven: vueApp.appData.uniqueAxisTransformationEven,
                uniqueAxisTransformationOdd: vueApp.appData.uniqueAxisTransformationOdd,
                layerMassAmu: vueApp.appData.layerMassAmu,
                layerSurfaceAng2: vueApp.appData.layerSurfaceAng2,
                numLayersBulk: vueApp.appData.numLayersBulk,
                maxLayers: vueApp.maxLayers
            };

            this.currentRequest = $.post({
                url: "/compute/api/modes/",
                //method: "POST",
                dataType: 'json',
                data: JSON.stringify(dataToSend)
            }).done(res => {
                vueApp.updateData(res);
			}).fail(function (jqXHR, textStatus, errorThrown) {
                if (textStatus === 'abort') {
                    // This is an aborted request
                    // console.log('ABORT', jqXHR);
                    return;
                }
                var warningMessage = 'Failed fetching the data: ' + textStatus;
                if (errorThrown) {
                    warningMessage += " (" + errorThrown + ": " + jqXHR.responseText + ")"; 
                }
                vueApp.warningMessage = warningMessage;
                vueApp.clearData();
            }).always( res => {
                vueApp.showSpinner = false;
                vueApp.currentRequest = null;
            });
        },
        updateData: function(data) {

            var vueApp = this;

            // NOTE! Data arrives in meV from the backend and is internally stored in meV in seriesWithData
            // Units conversion is done only while plotting (for the frequency units)
            vueApp.seriesWithData = _.chain([
                data.x,
                data.y_meV,
                data.isBackScattering,
                data.isRamanActive,
                data.isInfraredActive,
                data.irrepNames,
                data.alongX,
                data.alongY,
                data.alongZ])
            .unzip()
            .groupBy(function(point) {
                return (
                    "IR_" + (point[4] ? "Y" : "N") + 
                    "+Raman_" + (point[3] ? "Y" : "N") + 
                    "+BackScattering_" + (point[2] ? "Y" : "N")); 
            })
            .mapObject(function (val, key) {
                // val is of the form [
                //    [x1,y1,backscattering1, raman1, infrared1],
                //    [x2,y2,backscattering2, raman2, infrared2],
                //    ...
                // ]
                return _.map(val, function(singlePoint) {
                    return {x: singlePoint[0], y: singlePoint[1], 
                        meta: {
                            isBackScattering: singlePoint[2],
                            isRamanActive: singlePoint[3],
                            isInfraredActive: singlePoint[4],
                            irrepName: singlePoint[5],
                            alongX: singlePoint[6],
                            alongY: singlePoint[7],
                            alongZ: singlePoint[8]
                        }};
                })
            })
            .value();

            // Call redrawData, pass the value
            vueApp.redrawData();
        },
        onChangeFilter: function(event) {
            // The dropdown was changed. Update the internal variable, then call redrawData
            var vueApp = this;
            if (event.target.id == "modes-filter-x") {
                vueApp.modesFilterX = event.target.checked;
            } else if (event.target.id == "modes-filter-y") {
                vueApp.modesFilterY = event.target.checked;
            } else if (event.target.id == "modes-filter-z") {
                vueApp.modesFilterZ = event.target.checked;
            }
            // ignore other cases - it should not happen
            vueApp.redrawData();
        },
        onChangeUnits: function(event) {
            var vueApp = this;
            if (event.target.id == "frequency-units") {
                vueApp.frequencyUnits = event.target.value;
            }
            vueApp.redrawData();
        },
        filterData: function(data, filterX, filterY, filterZ) {
            // I am getting both the data and the modesFilter (the latter would also be in `this`,
            // but I defined this function as a "class method" instead)

            // Show only data points with non-zero component along a given direction
            return _.filter(data, function(point) {
                return (
                    (filterX && point.meta.alongX ) ||
                    (filterY && point.meta.alongY ) ||
                    (filterZ && point.meta.alongZ )
                );
            });
        },
        redrawData: function() {            
            var vueApp = this;

            // This is useful because it also keeps the order of the series
            var seriesNames = _.pluck(vueApp.seriesMetadata, 'internalName');

            var unitsInfo = energyUnitsMapping[vueApp.frequencyUnits];
            if (typeof(unitsInfo) == 'undefined') {
                unitsInfo = {
                    label: 'a.u.',
                    from_meV_factor: 1.,
                }
            }

            // clean-up all series
            while(vueApp.chart.series.length > 0) {
                vueApp.chart.series[0].remove();
            }
            
            var yPrecision = 4;  // Precision for y-value (the value itself, but mostly for the tooltip)

            var seriesWithDataNewUnits = _.mapObject(
                vueApp.seriesWithData, function(series, seriesName) {
                    // return the new value, for the same seriesName
                    return _.map(series, function(dataPoint) {
                        var newDataPoint = {};
                        Object.assign(newDataPoint, dataPoint);
                        // the original data point is assumed to be always returned in meV by the API
                        newDataPoint.y = newDataPoint.y * unitsInfo.from_meV_factor;
                        // Round the y value
                        newDataPoint.y = Math.round(newDataPoint.y * Math.pow(10, yPrecision)) / Math.pow(10, yPrecision);
                        return newDataPoint;
                    });
                }
            );

            // Fix the max of the y axis independent of the filter
            var maxY = _.max( // get the maxY from the maxYPerEachSeries
                _.map(
                    _.values(seriesWithDataNewUnits), // for each series...
                    function(series) {
                        //  ...get the max Y for this series
                        return _.max(series, function(dataPoint) {return dataPoint.y;}).y;
                    }
                )
            );

            // Recreate each of the series, in the right order
            _.each(
                seriesNames, function(name, seriesIdx) {
                    // Hide first, if it's going to be invisible anyway
                    if (_.has(seriesWithDataNewUnits, name)) {
                        var metaSeriesIndex = _.findIndex(vueApp.seriesMetadata, function(series){ return series.internalName == name; });

                        if (metaSeriesIndex == -1) {
                            console.warn("Unable to find series " + name);
                        }
                        else {
                            vueApp.chart.addSeries(
                                _.extend(
                                    vueApp.seriesMetadata[metaSeriesIndex], 
                                    {
                                        data: vueApp.filterData(seriesWithDataNewUnits[name], vueApp.modesFilterX, vueApp.modesFilterY, vueApp.modesFilterZ),
                                        showInLegend: true
                                    }
                                )
                            );
                        }
                    }
                }
            );

            // Fix the y range so it does not depend on the filtering
            vueApp.chart.yAxis[0].update({min: 0, max: maxY}); // The min is always zero
            vueApp.chart.yAxis[0].setTitle({text: "Frequency ω [" + unitsInfo.label + "]"});
        },
        clearData: function() {
            while(this.chart.series.length > 0) {
                this.chart.series[0].remove();
            }
        }, 
        createChart: function () {
            var _this = this;
            return Highcharts.chart('plotContainer', {
                chart: {
                    type: 'scatter',
                    zoomType: 'xy',
                    animation: false,
                    events: {
                        redraw: function (event) {
                            var xmin = event.target.axes[0].min;
                            var xmax = event.target.axes[0].max;
                        }
                    }
                },
                drilldown: {
                    animation: {
                        duration: 0
                    }
                },
                credits: {
                    enabled: false
                },        
                title: {
                    text: 'Modes' //, useHTML: true
                },
                /* subtitle: {
                    text: 'Move the slider to change the value of <emph>k</emph>',
                    useHTML: true

                }, */
                xAxis: {
                    title: {
                        text: 'Number of layers'
                    }
                },
                yAxis: {
                    title: {
                        text: 'Frequency ω'
                    }
                },
                tooltip: {
                    crosshairs: true,
                    shared: true,
                    useHTML: true
                },
                plotOptions: {
                    series: {
                        animation: {
                            duration: 0
                        }
                    },
                    scatter: {
                        marker: {
                            radius: 5,
                            states: {
                                hover: {
                                    enabled: true,
                                    lineColor: 'rgb(100,100,100)'
                                }
                            }
                        },
                        states: {
                            hover: {
                                marker: {
                                    enabled: false
                                }
                            }
                        },
                        tooltip: {
                            headerFormat: '<b>{series.name}</b><br>',
                            pointFormat: 'N={point.x}<br>Frequency={point.y}<br>Irrep name: {point.meta.irrepName}'
                        }
                    }
                },
                legend: {
                    /* layout: 'vertical',
                    align: 'left',
                    verticalAlign: 'top',
                    x: 100,
                    y: 70, 
                    floating: true, */
                    align: 'center',
                    verticalAlign: 'bottom',
                    x: 0,
                    y: 0,
                    backgroundColor: (Highcharts.theme && Highcharts.theme.legendBackgroundColor) || '#FFFFFF',
                    borderWidth: 1
                },
                series: [] // These will be created later
                    /* _.map(this.seriesMetadata, function (series) {
                        return _.omit(series, 'internalName'); // remove internalName 
                    }) */
            });
        }        
    },
    mounted: function () {
        // define immediately the debounced function (even before having the DOM element in-document)
        this.updateGraphDebounced = this.getUpdateDebouncedFunction();

        this.$nextTick(function () {
            // code that assumes this.$el is in-document

            // at this point everything should be loaded
            this.latexCode = window.appData.forceConstants.description;
            var vueApp = this;
            _.each(
                window.appData.forceConstants.variables, function(variable) {
                    vueApp.forceConstantVariables.push(variable);        
                }
            );

            // Define custom symbols
            window.Highcharts.SVGRenderer.prototype.symbols.mycrosscircle = function (x, y, w, h) {
                // I do not want to make a cross in a square "w x h", but
                // I want it to fit inside the circle inscribed in this rectangle/square.
                var squeeze = (Math.sqrt(2) - 1)/(2. * Math.sqrt(2)); 

                return [
                    'M', x + w * squeeze, y + h * squeeze,
                    'L', x + w * (1 - squeeze), y + h * (1 - squeeze),
                    'M', x + w * (1 - squeeze), y + h * squeeze,
                    'L', x + w * squeeze, y + h * (1 - squeeze),
                    'z',
                    'M', x, y+h/2, 'a', w/2, h/2, 0, 1, 0, w, 0, 'a', w/2, h/2, 0, 1, 0, -w, 0, 'z'
                ];
            };
            if (window.Highcharts.VMLRenderer) {
                window.Highcharts.VMLRenderer.prototype.symbols.mycrosscircle = window.Highcharts.SVGRenderer.prototype.symbols.mycrosscircle;
            }


            window.Highcharts.SVGRenderer.prototype.symbols.cross = function (x, y, w, h) {
                // I do not want to make a cross in a square "w x h", but
                // I want it to fit inside the circle inscribed in this rectangle/square.
                var squeeze = (Math.sqrt(2) - 1)/(2. * Math.sqrt(2)); 

                return [
                    'M', x + w * squeeze, y + h * squeeze,
                    'L', x + w * (1 - squeeze), y + h * (1 - squeeze),
                    'M', x + w * (1 - squeeze), y + h * squeeze,
                    'L', x + w * squeeze, y + h * (1 - squeeze),
                    'z'];
            };
            if (window.Highcharts.VMLRenderer) {
                window.Highcharts.VMLRenderer.prototype.symbols.cross = window.Highcharts.SVGRenderer.prototype.symbols.cross;
            }
            window.Highcharts.SVGRenderer.prototype.symbols.mycircle = function (x, y, w, h) {
               /* d="
                M (CX - R), CY
                a R,R 0 1,0 (R * 2),0
                a R,R 0 1,0 -(R * 2),0
                z
                "*/
                return ['M', x, y+h/2, 'a', w/2, h/2, 0, 1, 0, w, 0, 'a', w/2, h/2, 0, 1, 0, -w, 0, 'z'];
            };
            if (window.Highcharts.VMLRenderer) {
                window.Highcharts.VMLRenderer.prototype.symbols.mycircle = window.Highcharts.SVGRenderer.prototype.symbols.mycircle;
            }

            // Create the chart, store in app
            this.chart = this.createChart();

            window.chart = this.chart;

            // Fetch data immediately, without debouncing
            this.fetchAndUpdateData(); 
            }
        )
      }
})



