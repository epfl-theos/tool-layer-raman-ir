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
        appData: appData
    },
    methods: {
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
            this.currentRequest = $.post({
                url: "http://localhost:5300/api/modes/",
                //method: "POST",
                dataType: 'json',
                data: JSON.stringify({
                    forceConstantParams: _.object(
                        _.map(vueApp.appData.forceConstants.variables, function(variable) {return variable['name'];}),
                        _.map(vueApp.appData.forceConstants.variables, function(variable) {return variable['value'];})
                    ),
                    matrices: vueApp.appData.forceConstants.matrices,
                    symmetryInfo: vueApp.appData.symmetryInfo,
                    maxLayers: vueApp.maxLayers,
                })
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
            console.log(data);
            this.chart.series[0].setData(_.zip(data.x, data.y)); // sin
            // this.chart.setTitle({text: 'Modes (<emph>k</emph> = ' + data.k + ')', useHTML: true});
        },
        clearData: function() {
            this.chart.series[0].setData([[]]); 
        },
        createChart: function () {
            return Highcharts.chart('plotContainer', {
                chart: {
                    type: 'scatter',
                    zoomType: 'xy',
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
                            pointFormat: 'x={point.x}, y={point.y}'
                        }
                    }
                },
                legend: {
                    layout: 'vertical',
                    align: 'left',
                    verticalAlign: 'top',
                    x: 100,
                    y: 70,
                    floating: true,
                    backgroundColor: (Highcharts.theme && Highcharts.theme.legendBackgroundColor) || '#FFFFFF',
                    borderWidth: 1
                },
                series: [{
                    name: 'Modes',
                    marker: {
                        symbol: 'square'
                    },
                    color: 'rgba(223, 83, 83, .5)',
                    data: [],
                    dataLabels: {
                        useHTML: true
                    }
                }, {
                    name: 'Second series',
                    color: 'rgba(119, 152, 191, .5)',
                    marker: {
                        symbol: 'diamond'
                    },
                    data: []
                }]
            });
        }        
    },
    mounted: function () {
        // define immediately the debounced function (even before having the DOM element in-document)
        this.updateGraphDebounced = this.getUpdateDebouncedFunction();

        this.$nextTick(function () {
            // code that assumes this.$el is in-document
            //console.log('mounted', this.$el);

            // Create the chart, store in app
            this.chart = this.createChart();

            // Fetch data immediately, without debouncing
            this.fetchAndUpdateData(); 
            }
        )
      }
})



