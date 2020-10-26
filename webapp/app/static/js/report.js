var dataObject = null;
Smits.PhyloCanvas.Render.Parameters.Rectangular["paddingY"] = 40;
Smits.PhyloCanvas.Render.Parameters.Rectangular["paddingX"] = 90;
function renderTree(dataObject,width,height) {
$('#svgCanvas').empty();
$('#sizebtns').attr("style","display:block");
phylocanvas = new Smits.PhyloCanvas(
				dataObject,
				'svgCanvas',
				width, height
				);
				/* var canvasHeight = parseInt($('svg').attr("height"));
				var canvasWidth = parseInt($('svg').attr("width")); */
				$('tspan:contains(QS--)').attr("fill","#0000ff");
				$('tspan:contains(TS--)').attr("fill","#1f7a1f");
				$('tspan:contains(OG--)').attr("fill", '#cc3300');
/*				$('tspan').each(function(objindex,obj) {
				console.log($(obj));
				console.log(Object.keys($(obj)));
				if ($(obj).html().substring(0,3) == "___") {
				$(obj).attr("fill","#0000ff");
				}
				}); */
                var treeSvg = $('#svgCanvas').html();
                var fileNameToSave = "treetest3.svg";
                treeSvg = '<?xml version="1.0" standalone="no"?>\r\n' + treeSvg;
                var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(treeSvg);
                $('#imgdownl').attr("href",url);
                $('#imgdownl').attr("download", fileNameToSave);
}
function repairSize() {
var canvasHeight = parseInt($('svg').attr("height"));
var canvasWidth = parseInt($('svg').attr("width"));
console.log(canvasWidth);
$('text:has(tspan)').each(function(objindex, obj) {
    var textHeight = parseInt($(this).attr("y"));
    var textX = parseInt($(this).attr("x"));
    var textLength = $(this).children("tspan").text().length;
    var textWidth = textX + (7*textLength);
    if (textHeight > canvasHeight) {
        console.log(textHeight);
        $('svg').attr("height",textHeight +10);
        $('#svgCanvas').attr("height",textHeight + 10);
        }
    if (textWidth > canvasWidth) {
        $('svg').attr("width",textWidth + 30);
        $('#svgCanvas').attr("width",textWidth + 30);
    }

});
}

function resizeTree(resizeDir) {
var canvasHeight = parseInt($('svg').attr("height"));
var canvasWidth = parseInt($('svg').attr("width"));
if (resizeDir == "upV") {
    renderTree(dataObject,canvasWidth,canvasHeight + 100);
    /*$('#svgCanvas').attr("height", canvasHeight+100 + Math.round((2*((canvasHeight+100)/100))));
    $('svg').attr("height", canvasHeight+100 + Math.round((2*((canvasHeight+100)/100))));*/
    } else if (resizeDir == "downV") {
    renderTree(dataObject,canvasWidth,canvasHeight - 100);
    /*$('#svgCanvas').attr("height", canvasHeight - 100 + Math.round((2*((canvasHeight - 100)/100))));
    $('svg').attr("height", canvasHeight - 100 + Math.round((2*((canvasHeight - 100)/100)))); */
    } else if (resizeDir == "upH") {
    renderTree(dataObject,canvasWidth + 100, canvasHeight);
    } else if (resizeDir == "downH") {
    renderTree(dataObject,canvasWidth - 100, canvasHeight);
    }
    $('text:has(tspan)').each(function(objindex, obj) {
        var textHeight = parseInt($(this).attr("y"));
        var textX = parseInt($(this).attr("x"));
        var textLength = $(this).children("tspan").text().length;
        var textWidth = textX + (7*textLength);
        if (textWidth > canvasWidth) {
            console.log(textWidth);
            repairSize();
            }
        });
}

function treeSuccess(data,textStatus,xhr) {
if (data != "false") {
dataObject = {
				newick: data,
				fileSource: true
			};
			renderTree(dataObject,1000,1000);
			repairSize();
			clearInterval(timer);
}
}

function treeError(xhr,ajaxOptions,thrownError) { // communication error
    console.log(xhr,ajaxOptions,thrownError);
}
function drawTree() {
var jobid = $('#jobinfo').val();
//console.log(jobid);
$.ajax({
        url: '/results/'+jobid+'/tree',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: treeSuccess,
        error: treeError});
}

$("button[data-toggle='tooltip']").on('click',function(){
    $(this).blur();
})

var timer = setInterval(drawTree, 3000);