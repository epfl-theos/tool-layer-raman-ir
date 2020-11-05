function toggleStrVisInteraction(enableStrInteraction, str_overlay_id){
    if (enableStrInteraction){
        // enable interaction here
        $("#"+str_overlay_id).css("display", "none")
        .css("-webkit-touch-callout", "auto")
        .css("-webkit-user-select", "auto")
        .css("-khtml-user-select", "auto")
        .css("-moz-user-select", "auto")
        .css("-ms-user-select", "auto")
        .css("user-select", "auto");
    }
    else{
        // disable interaction here
        $("#"+str_overlay_id).css("display", "table")
        .css("-webkit-touch-callout", "none")
        .css("-webkit-user-select", "none")
        .css("-khtml-user-select", "none")
        .css("-moz-user-select", "none")
        .css("-ms-user-select", "none")
        .css("user-select", "none");
    }
};


function postLoadByApplet(viewerApplet){
    // This function should be called every time a `load` command is issued, to recreate the correct bonds
    Jmol.script(viewerApplet, "select *; connect delete; " + jsmolBondCommand);
    // Hide the bonds based on the state of the checkbox
    if (viewerApplet._id.startsWith('jmolApplet')) {
        showBonds(viewerApplet._id.substr('jmolApplet'.length));
    }
};

function postLoadByNameSuffix(viewerNameSuffix){
    // Utility function to call postLoadByApplet when only the suffix of the viewer name is known
    postLoadByApplet(eval('jmolApplet' + viewerNameSuffix));
};


function jsmolCrystal(data, parentHtmlId, appletName, supercellOptions) {
    var parentDiv = document.getElementById(parentHtmlId);
    var the_width = parentDiv.offsetWidth - 5;
    var the_height = parentDiv.offsetHeight - 5;

    var Info = {
        width: the_width,
        height: the_height,
        debug: false,
        color: "#FFFFFF",
        use: "HTML5",
        j2sPath: "../../user_static/js/jsmol/j2s",
        serverURL: "../../user_static/js/jsmol/php/jsmol.php",
        console: appletName + "_infodiv"
    };

    var jsmolStructureviewer = Jmol.getApplet(appletName, Info);

    // I `set autobond off` so bonds are not computed; remember, after *any* load command, to re-call the postLoadByApplet function
    // (or the wrapper function postLoadByNameSuffix)
    if (supercellOptions === undefined){
        var loadingScript = 'color cpk; set autobond off; load INLINE "' + data + '"; wireframe 0.15; spacefill 23%';
    } else {
        var loadingScript = 'color cpk; set autobond off; load INLINE "' + data + '" {' + supercellOptions[0] + ' ' + supercellOptions[1] + ' ' + supercellOptions[2] + '}; wireframe 0.15; spacefill 23%';
    }

    //draw x,y,z axes instead of a,b,c vectors as default
    loadingScript+= '; axes off; draw xaxis ">X" vector {0 0 0} {2 0 0} color red width 0.15; draw yaxis ">Y" vector {0 0 0} {0 2 0} color green width 0.15; draw zaxis ">Z" vector {0 0 0} {0 0 2} color blue width 0.15';

    //do not show info on top left
    loadingScript+= '; unitcell primitive';

    //Sets the unit cell line diameter in Angstroms
    loadingScript+= '; unitcell 2';

    // antialiasDisplay ON
    loadingScript+= '; set antialiasDisplay on';

    //Zooms to the setting that fills the screen with the currently displayed atoms
    loadingScript+= "; set zoomLarge false";

    //remove JSmol logo
    loadingScript+= '; set frank off';

    Jmol.script(jsmolStructureviewer, loadingScript);

    postLoadByApplet(jsmolStructureviewer);

    //parentDiv.innerHTML = Jmol.getAppletHtml(jsmolStructureviewer);
    //$("#" + parentHtmlId).html(Jmol.getAppletHtml(jsmolStructureviewer));

    return jsmolStructureviewer;
}

var cellLine = "; unitcell 2";

function toggleRotation(viewerNameSuffix) {
    if ($("#spin-input-" + viewerNameSuffix).is(":checked")){
        var jmolscript = "spin on";
    } else {
        var jmolscript = "spin off";
    }
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);
    return jmolscript;
};

function showBonds(viewerNameSuffix) {
    if ($("#bonds-input-" + viewerNameSuffix).is(":checked")){
        var jmolscript = "wireframe 0.15";
    } else {
        var jmolscript = "wireframe off";
    }
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);
    return jmolscript;
};

function showPacked(viewerNameSuffix) {
    var nx = $('#nx-' + viewerNameSuffix).val();
    var ny = $('#ny-' + viewerNameSuffix).val();
    var nz = $('#nz-' + viewerNameSuffix).val();

    if ($("#packed-input-" + viewerNameSuffix).is(":checked")){
        var jmolscript = "save orientation 0; load '' {" + nx + " " + ny + " " + nz + "} packed; unitcell primitive; restore orientation 0" + jsmolDrawAxes(viewerNameSuffix) + cellLine + "; " + showLabels(viewerNameSuffix) + "; " + showBonds(viewerNameSuffix);
    } else {
        var jmolscript = "save orientation 0; load '' {" + nx + " " + ny + " " + nz + "}; unitcell primitive; restore orientation 0" + jsmolDrawAxes(viewerNameSuffix) + cellLine + "; " + showLabels(viewerNameSuffix) + "; " + showBonds(viewerNameSuffix);
    }
    $("#spin-input-" + viewerNameSuffix).prop("checked", false);
    $("#spheres-input-" + viewerNameSuffix).prop("checked", false);
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);

    postLoadByNameSuffix(viewerNameSuffix);

    return jmolscript;
};

function showLabels(viewerNameSuffix) {
    if ($("#labels-input-" + viewerNameSuffix).is(":checked")){
        var jmolscript = "label %a";
    } else {
        var jmolscript = "label off";
    }
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);
    return jmolscript;
};

function showSpheres(viewerNameSuffix) {
    if ($("#spheres-input-" + viewerNameSuffix).is(":checked")){
        var jmolscript = "spacefill on";
    } else {
        var jmolscript = "spacefill 23%";
    }
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);
    return jmolscript;
};

function jsmolDrawAxes(viewerNameSuffix) {
    var e = document.getElementById("axesMenu-" + viewerNameSuffix);
    var selectedAxes = e.options[e.selectedIndex].value;
    switch (selectedAxes){
        case "xyz":
            var jmolscript = "; axes off; draw xaxis '>X' vector {0 0 0} {2 0 0} color red width 0.15; draw yaxis '>Y' vector {0 0 0} {0 2 0} color green width 0.15; draw zaxis '>Z' vector {0 0 0} {0 0 2} color blue width 0.15";
            break;
        case "abc":
            var jmolscript = "; draw xaxis delete; draw yaxis delete; draw zaxis delete; set axesMode 2; axes 5";
            break;
        case "noaxes":
            var jmolscript = "; draw xaxis delete; draw yaxis delete; draw zaxis delete; axes off";
    }
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);
    return jmolscript;
};

function jsmolSupercell(viewerNameSuffix) {
    var nx = $('#nx-' + viewerNameSuffix).val();
    var ny = $('#ny-' + viewerNameSuffix).val();
    var nz = $('#nz-' + viewerNameSuffix).val();
    $("#spin-input-" + viewerNameSuffix).prop("checked", false);
    $("#spheres-input-" + viewerNameSuffix).prop("checked", false);
    $("#packed-input-" + viewerNameSuffix).prop("checked", false);
    var jmolscript = "save orientation 0; load '' {" + nx + " " + ny + " " + nz + "}; unitcell primitive; restore orientation 0" + jsmolDrawAxes(viewerNameSuffix) + cellLine + "; " + showLabels(viewerNameSuffix) + "; " + showBonds(viewerNameSuffix);
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), jmolscript);

    postLoadByNameSuffix(viewerNameSuffix);
};

function jsmolResetSupercell(viewerNameSuffix, nx, ny, nz) {
    $("#spin-input-" + viewerNameSuffix).prop("checked", false);
    $("#spheres-input-" + viewerNameSuffix).prop("checked", false);
    $("#packed-input-" + viewerNameSuffix).prop("checked", false);
    // reset nx, ny, nz to default values
    $('#nx-' + viewerNameSuffix).val(nx);
    $('#ny-' + viewerNameSuffix).val(ny);
    $('#nz-' + viewerNameSuffix).val(nz);
    Jmol.script(eval('jmolApplet' + viewerNameSuffix),
    "save orientation 0; load '' {" + nx + " " + ny + " " + nz + "}; unitcell primitive; restore orientation 0" + jsmolDrawAxes(viewerNameSuffix) + cellLine + "; " + showLabels(viewerNameSuffix) + "; " + showBonds(viewerNameSuffix));

    postLoadByNameSuffix(viewerNameSuffix);
};

function centerXaxis(viewerNameSuffix){
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), "moveto 1 axis x");
};

function centerYaxis(viewerNameSuffix){
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), "moveto 1 axis y");
};

function centerZaxis(viewerNameSuffix){
    Jmol.script(eval('jmolApplet' + viewerNameSuffix), "moveto 1 axis z");
};


$.fn.bindFirst = function(name, fn) {
    var elem, handlers, i, _len;
    this.bind(name, fn);
    for (i = 0, _len = this.length; i < _len; i++) {
      elem = this[i];
      handlers = jQuery._data(elem).events[name.split('.')[0]];
      handlers.unshift(handlers.pop());
    }
  };

function enableDoubleTap(element, callback, ignoreOnMove) {
    /* Enable double-tap event for phones */
    element.dbltapTimeout = undefined;
    element.shortTap = false;

    var preventOnMove = false;
    if (typeof ignoreOnMove !== 'undefined') {
        preventOnMove = ignoreOnMove;
    }

    // Manual detect of double tap
    //element.addEventListener('touchend', function(event) {
    $(element).bindFirst('touchend', function(event) {
      if (typeof element.dbltapTimeout !== 'undefined') {
          // start disabling any timeout that would reset shortTap to false
          clearTimeout(element.dbltapTimeout);
      }
      if (element.shortTap) {
          // if here, there's been another tap a few ms before
          // reset the variable and do the custom action
          element.shortTap = false;
          event.preventDefault();
          event.stopImmediatePropagation();
          callback();
      }
      else {
          if (event.targetTouches.length != 0) {
              // activate this only when there is only a finger
              // if more than one finger is detected, cancel detection
              // of double tap
              if (typeof element.dbltapTimeout !== 'undefined') {
                  // disable the timeout
                  clearTimeout(element.dbltapTimeout);
                  element.shortTap = false;
              }
            return;
          }
          // If we are here, no tap was recently detected
          // mark that a tap just happened, and start a timeout
          // to reset this
          element.shortTap = true;

          element.dbltapTimeout = setTimeout(function() {
              // after 500ms, reset shortTap to false
              element.shortTap = false;
          }, 500);
      }
  });
  element.addEventListener('touchcancel', function(event) {
      if (typeof element.dbltapTimeout !== 'undefined') {
          // disable the timeout if the touch was canceled
          clearTimeout(element.dbltapTimeout);
          element.shortTap = false;
      }
  });
  if (!preventOnMove) {
        element.addEventListener('touchmove', function(event) {
        if (typeof element.dbltapTimeout !== 'undefined') {
            // disable the timeout if the finger is being moved
            clearTimeout(element.dbltapTimeout);
            element.shortTap = false;
        }
    });
  }
}
