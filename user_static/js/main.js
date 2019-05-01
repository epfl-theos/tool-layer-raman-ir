'use strict';

// Get the data from the static information compiled into the HTML page from the server;
// Also, prevent these properties to be changed 
// Object.freeze(appData); 

var app = new Vue({
    el: '#app',
    data: {
        maxLayers: 10,
        warningMessage: '',
        showSpinner: false,
        currentRequest: null,
        appData: appData,
        latexCode: "",
        forceConstantVariables: [],

        seriesMetadata: Object.freeze([ // Initial series definition
            {
                internalName: 'IR_Y+Raman_Y+BackScattering_Y',
                name: 'Infrared + Raman active (visible in backscattering)',
                marker: {
                    symbol: 'square'
                },
                color: 'rgba(223, 83, 83, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_Y+Raman_Y+BackScattering_N',
                name: 'Infrared + Raman active (not visible in backscattering)',
                marker: {
                    symbol: 'square'
                },
                color: 'rgba(0, 0, 0, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_Y+Raman_N+BackScattering_Y',
                name: 'Infrared active (visible in backscattering)',
                marker: {
                    symbol: 'triangle'
                },
                color: 'rgba(223, 83, 83, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_Y+Raman_N+BackScattering_N',
                name: 'Infrared active (not visible in backscattering)',
                marker: {
                    symbol: 'triangle'
                },
                color: 'rgba(0, 0, 0, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_Y+BackScattering_Y',
                name: 'Raman active (visible in backscattering)',
                marker: {
                    symbol: 'circle'
                },
                color: 'rgba(223, 83, 83, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_Y+BackScattering_N',
                name: 'Raman active (not visible in backscattering)',
                marker: {
                    symbol: 'circle'
                },
                color: 'rgba(0, 0, 0, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_N+BackScattering_Y',
                name: 'Inactive (visible in backscattering)',
                marker: {
                    symbol: 'square'
                },
                color: 'rgba(223, 83, 83, .5)',
                data: [],
                dataLabels: {
                    useHTML: true
                }    
            },
            {
                internalName: 'IR_N+Raman_N+BackScattering_N',
                name: 'Inactive (not visible in backscattering)',
                marker: {
                    symbol: 'square'
                },
                color: 'rgba(0, 0, 0, .5)',
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
                symmetryInfo: vueApp.appData.symmetryInfo,
                maxLayers: vueApp.maxLayers,
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
                    warningMessage += " (" + errorThrown + ")";  
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

            // This is useful because it also keeps the order of the series
            var seriesNames = _.pluck(vueApp.seriesMetadata, 'internalName');

            // will eventually contain all data, for now put empty defaults
            // form: 
            /*
                {
                    'IR_Y+Raman_Y+BackScattering_Y': {
                        data: [[]],
                        visible: false
                    },
                    'IR_Y+Raman_Y+BackScattering_N': {
                        data: [[]],
                        visible: false
                    },
                    ...
                }
            */

            // construct a dictionary with *only* the data series you got from the API
            // form: 
            /*
                {
                    'IR_Y+Raman_Y+BackScattering_Y': [[x1,y1], [x2,y2]],
                    'IR_Y+Raman_Y+BackScattering_N': [[x1,y1], [x2,y2]],
                    ...
                }
            */
            var seriesWithData = _.chain([data.x, data.y, data.isBackScattering, data.isRamanActive, data.isInfraredActive])
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
                            isInfraredActive: singlePoint[4]
                        }};
                })
            })
            .value();
            
            //console.log(seriesNames, seriesWithData);

            // clean-up all series
            while(vueApp.chart.series.length > 0) {
                vueApp.chart.series[0].remove();
            }

            // Recreate each of the series, in the right order
            _.each(
                seriesNames, function(name, seriesIdx) {
                    // Hide first, if it's going to be invisible anyway
                    if (_.has(seriesWithData, name)) {
                        var metaSeriesIndex = _.findIndex(vueApp.seriesMetadata, function(series){ return series.internalName == name; });

                        if (metaSeriesIndex == -1) {
                            console.warn("Unable to find series " + name);
                        }
                        else {
                            vueApp.chart.addSeries(
                                _.extend(
                                    vueApp.seriesMetadata[metaSeriesIndex], 
                                    {
                                        data: seriesWithData[name],
                                        showInLegend: true
                                    }
                                )
                            );
                        }
                    }
                }
            );

            // this.chart.setTitle({text: 'Modes (<emph>k</emph> = ' + data.k + ')', useHTML: true});
        },
        clearData: function() {
            while(this.chart.series.length > 0) {
                this.chart.series[0].remove();
            }
        },
        createChart: function () {
            return Highcharts.chart('plotContainer', {
                chart: {
                    type: 'scatter',
                    zoomType: 'xy'
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
                        text: 'Frequency ω [a.u.]'
                    }
                },
                tooltip: {
                    crosshairs: true,
                    shared: true,
                    useHTML: true
                },
                plotOptions: {
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
                            pointFormat: 'N={point.x}<br>Frequency={point.y}<br>isBackScattering={point.meta.isBackScattering}'
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
            //console.log('mounted', this.$el);

            // at this point everything should be loaded
            this.latexCode = window.appData.forceConstants.description;
            var vueApp = this;
            _.each(
                window.appData.forceConstants.variables, function(variable) {
                    vueApp.forceConstantVariables.push(variable);        
                }
            );

            // Define custom symbols
            window.Highcharts.SVGRenderer.prototype.symbols.cross = function (x, y, w, h) {
                //console.log('cross', x, y, w, h);
                return ['M', x, y, 'L', x + w, y + h, 'M', x + w, y, 'L', x, y + h, 'z'];
            };
            if (window.Highcharts.VMLRenderer) {
                window.Highcharts.VMLRenderer.prototype.symbols.cross = window.Highcharts.SVGRenderer.prototype.symbols.cross;
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



