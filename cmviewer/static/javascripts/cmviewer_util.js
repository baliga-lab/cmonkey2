// These function were pulled out of index.html to make it more maintainable, need to clean
// up here when it's ready
var currentIteration = 1;

function updateIterationSelector() {
    $.ajax({ url: '/iteration_select', success: function(html) {
                 $(html.trim()).replaceAll('#iteration_select');
                 $('#select_iteration').change(
                     function (event) {
                         currentIteration = $(this).val();
                         var iteration = $(this).val();
                         // Update cluster list and networks
                         var activeTab = $('#tabs').tabs("option", "active");
                         $('#cluster-list').DataTable().ajax.url('/clusters/' + iteration).load();
                         clearClusterDetailView();

                         if (activeTab == 2) { // only update this if we are in the network tab
                             var residuals = $('#residual-slider').slider("values");
                             var evalues = $('#evalue-slider').slider("values");
                             initCytoweb(iteration, residuals[0], residuals[1], evalues[0], evalues[1]);
                         }
                     });
             }
           });
}

function updateRunStatus(pgbarSelector) {
    $.ajax({ url: '/run_status', success: function(data) {
                 var progress = parseFloat(data.progress);
                 $('#progressbar').progressbar("value", progress);
                 $('.progress-label').text(progress + "%");
             }
           });
    updateIterationSelector();
}

function assignClusterClickHandlers() {
    var cluster = $(this).attr('id');
    $.get('/cluster/' + currentIteration + '/' + cluster,
          function(html) {
              // we need the trim() call in order for jQuery 1.9.x to
              // recognize the output as valid HTML
              $(html.trim()).replaceAll('#cluster-view');
              // notify Firegoose to rescan the page
              var ev = document.createEvent("Events");
              ev.initEvent("gaggleDataEvent", true, false);
              document.dispatchEvent(ev);
          });
    return false;
}

function clearClusterDetailView() {
    $("<span id=\"cluster-view\">Please select a cluster</span>").replaceAll('#cluster-view');
}

function showClusterDetails(cluster) {
    $.get('/cluster/' + currentIteration + '/' + cluster,
          function(html) {
              $(html.trim()).replaceAll('#cluster-view');
              $("#tabs").tabs( "option", "active", 1);
              // notify Firegoose to rescan the page
              var ev = document.createEvent("Events");
              ev.initEvent("gaggleDataEvent", true, false);
              document.dispatchEvent(ev);
          });
}

function initCytoweb(iteration, minResidual, maxResidual, minEvalue, maxEvalue) {
    var cy = $('#cy').cytoscape(
        {
            style: cytoscape.stylesheet()
                .selector('node').css(
                    {
                        'content': 'data(name)',
                        'text-valign': 'center',
                        'color': 'black',
                        'min-zoomed-font-size': 7
                    })
                .selector('edge').css(
                    {
                        'target-arrow-shape': 'none',
                        'curve-style': 'haystack'
                    })
                .selector('.clusters').css(
                    {
                        'color': 'green', 'background-color': 'red',
                        'width': 10, 'height': 10, 'font-size': '8'
                    })
                .selector('.genes').css(
                    {
                        'width': 5, 'height': 5, 'font-size': '3'
                    })
                .selector('.motifs').css(
                    {
                        'width': 5, 'height': 5, 'font-size': '3', 'shape': 'triangle',
                        'background-color': 'green'
                    })
                .selector(':selected').css(
                    {
                        'background-color': 'black',
                        'line-color': 'black'
                    })
                .selector('.faded').css(
                    {
                        'opacity': 0.25,
                        'text-opacity': 0
                    }),
            ready: function() {
                var cy = this;
                cy.startBatch();
                jQuery.getJSON('/cytoscape_nodes/' + iteration,
                               {
                                   min_residual: minResidual,
                                   max_residual: maxResidual,
                                   min_evalue: minEvalue,
                                   max_evalue: maxEvalue
                               },
                               function(nodes) {
                                   cy.add(nodes);
                                   jQuery.getJSON('/cytoscape_edges/' + iteration,
                                                  {
                                                      min_residual: minResidual,
                                                      max_residual: maxResidual,
                                                      min_evalue: minEvalue,
                                                      max_evalue: maxEvalue
                                                  },
                                                  function(edges) {
                                                      cy.add(edges);
                                                      cy.endBatch();
                                                      cy.layout({name: 'springy', animate: false, fit: true});
                                                      cy.fit();
                                                      cy.on('tap', 'node',
                                                            function() {
                                                                // open tab 1 with a filled in detail view
                                                                if (this.hasClass('clusters')) {
                                                                    showClusterDetails(this.id());
                                                                }
                                                            });
                                                  });
                               });
            },
            pixelRatio: 1.0,
            motionBlur: true,
            hideEdgesOnViewport: true,
            hideLabelsOnViewport: true,
            textureOnViewport: true
        });
}
