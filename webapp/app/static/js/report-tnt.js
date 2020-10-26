var jobid = $('#jobinfo').val();

//  ANI group coloring of labels
function nodeFill(node,aniGroupCutoff) {
    if (aniGroupCutoff == '97') {
    var nodeAni = node.property('ANI97');
    if (nodeAni in aniTable[0]) {
        return aniTable[0][nodeAni];
    } else {
        return "black";
    }
    } else if (aniGroupCutoff == '95') {
    var nodeAni = node.property('ANI95');
    if (nodeAni in aniTable[1]) {
        return aniTable[1][nodeAni];
    } else {
        return "black";
    }
    } else if (aniGroupCutoff == '90') {
    var nodeAni = node.property('ANI90');
    if (nodeAni in aniTable[2]) {
        return aniTable[2][nodeAni];
    } else {
        return "black";
    }
    }
    else if (aniGroupCutoff == '85') {
    var nodeAni = node.property('ANI85');
    if (nodeAni in aniTable[3]) {
        return aniTable[3][nodeAni];
    } else {
        return "black";
    }
    }
    else if (aniGroupCutoff == '80') {
    var nodeAni = node.property('ANI80');
    if (nodeAni in aniTable[4]) {
        return aniTable[4][nodeAni];
    } else {
        return "black";
    }
    }

}
//  coloring single ANI groups
function subtreeFill(node, searchId,aniCutoff) { //see nodeFill; how to get ANI cutoff?
if (aniCutoff == 'ANI97') {
if (node.property('ANI97') == searchId) {
    if (searchId in aniTable[0]) {
        return aniTable[0][searchId];
    } else {
        return "black";
    }
} else {
    return "black";
}} else if (aniCutoff == 'ANI95') {
if (node.property('ANI95') == searchId) {
    if (searchId in aniTable[1]) {
        return aniTable[1][searchId];
    } else {
        return "black";
    }
} else {
    return "black";
}
} else if (aniCutoff == 'ANI90') {
if (node.property('ANI90') == searchId) {
    if (searchId in aniTable[2]) {
        return aniTable[2][searchId];
    } else {
        return "black";
    }
} else {
    return "black";
}
} else if (aniCutoff == 'ANI85') {

if (node.property('ANI85') == searchId) {
    if (searchId in aniTable[3]) {
        return aniTable[3][searchId];
    } else {
        return "black";
    }
} else {
    return "black";
}
} else if (aniCutoff == 'ANI80') {
if (node.property('ANI80') == searchId) {
    if (searchId in aniTable[4]) {
        return aniTable[4][searchId];
    } else {
        return "black";
    }
} else {
    return "black";
}
} else {
    return "black";
}
}
//  coloring of nodes as type strains/query seqs/outgroups; also used for coloring labels
function strainTypeNodes(node) {
if (node.node_name().search("QS--") != -1) {
              return '#0000ff';
            }
            else if (node.node_name().search("OG--") != -1) { // marking as outgroup currently has higher priority than marking as TS
              return '#cc3300';
            }
            else if (node.node_name().search("TS--") != -1) {
              return '#1f7a1f';
            }
            return '#888888';
}

// setting up SVG; needs to be called after every visual change to the tree that should be reflected in downloaded image
function setSvg() {
$('svg').attr('xmlns',"http://www.w3.org/2000/svg");
$('svg').attr('xmlns:xlink', "http://www.w3.org/1999/xlink");
var treeSvg = $('#svgCanvas').html();
                var fileNameToSave = "treetest3.svg";
                treeSvg = '<?xml version="1.0" standalone="no"?>\r\n' + treeSvg;
                var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(treeSvg);
                //console.log(url);
                $('#imgdownl').attr("href",url);
                $('#imgdownl').attr("download", fileNameToSave);
                }
function sortAniList(aniDictToSort) {
// create array of arrays storing highest branch height (lowest y-offset) & ANI id
      var sortableAni = [];
      for (var z in aniDictToSort) {
      if (aniDictToSort[z].length > 1) {
      higherBranch = Math.min.apply(null,aniDictToSort[z]);
      sortableAni.push([higherBranch, z]);
      }
      }
      //sort by y-position
      sortableAni.sort(function(a,b){
        return a[0]-b[0];
      });
      return sortableAni;
}
// actual tree rendering
function treeSuccess(data,textStatus,xhr) {
if (data != "false") {
	$('#svgCanvas').empty();
	//$('#svgCanvas').addClass("hidden");
    $('#sizebtns').attr("style","display:block");
    $('#searchdiv').attr("style","display:block");
    var newick = data.toString();
    //console.log(newick);
    tree = tnt.tree();
     tree
     .data (tnt.tree.parse_newick(newick))
    .node_display(tree.node_display()
            .size(function(node) {
            if (node.is_leaf()) {
                return 4;
            } else {
                return 2;
            }
            })
            .fill(function(node) {
            if (node.is_leaf()) {
                return strainTypeNodes(node);
            }
            return '#888888';
            })
        )
        .label (tnt.tree.label.text()
            .fontsize(12)
            .height(30)
            .text (function (node) {
            return node.node_name();
        })
            .color("black")
        )
        .layout(tnt.tree.layout.vertical()
            .width(1000)
            .scale(false)
        );
          tree(document.getElementById("svgCanvas")); // tree needs to be rendered at this point to get y-positions of nodes later
//          clearInterval(timer);
        var aniDict = [];
        aniDict[0] = {};
        aniDict[1] = {};
        aniDict[2] = {};
        aniDict[3] = {};
        aniDict[4] = {};
        /*aniDict[95] = {};
        aniDict[96] = {};
        aniDict[97] = {};
        aniDict[98] = {};
        aniDict[99] = {};*/
        root = tree.root();
        leaves = root.get_all_leaves();
        // get all ANI ids and y-positions
        for (var leafNode in leaves) {
        var nodeName = leaves[leafNode].node_name();
        var gcfIndex = nodeName.search('GCF_');
        if (gcfIndex != -1) {
        var nodeGcf = nodeName.slice(gcfIndex,gcfIndex+13);
        var nodeAni = [aniMasterList[nodeGcf][0], aniMasterList[nodeGcf][1], aniMasterList[nodeGcf][2], aniMasterList[nodeGcf][3],aniMasterList[nodeGcf][4]];
        leaves[leafNode].property('ANI97', nodeAni[0]); // store as property so the name won't have to be searched every time
        leaves[leafNode].property('ANI95', nodeAni[1]);
        leaves[leafNode].property('ANI90', nodeAni[2]);
        leaves[leafNode].property('ANI85', nodeAni[3]);
        leaves[leafNode].property('ANI80', nodeAni[4]);
        var nodeHeight = $('text.tnt_tree_label:contains("'+ nodeName+'")').offset().top; //get y-position of branches relative to page
        for (var b in nodeAni) {
                if (!(nodeAni[b] in aniDict[b])) {
                    aniDict[b][nodeAni[b]] = [];
                }
                aniDict[b][nodeAni[b]].push(nodeHeight);
        }
        //long version:
        /*if (!(nodeAni[0] in aniDict[95])) {
            aniDict[95][nodeAni[0]] = [];
        }
        aniDict[95][nodeAni[0]].push(nodeHeight);
        if (!(nodeAni[1] in aniDict[96])) {
            aniDict[96][nodeAni[1]] = [];
        }
        aniDict[96][nodeAni[1]].push(nodeHeight);
        if (!(nodeAni[2] in aniDict[97])) {
            aniDict[97][nodeAni[2]] = [];
        }
        aniDict[97][nodeAni[2]].push(nodeHeight);
        if (!(nodeAni[3] in aniDict[98])) {
            aniDict[98][nodeAni[3]] = [];
        }
        aniDict[98][nodeAni[3]].push(nodeHeight);
        if (!(nodeAni[4] in aniDict[99])) {
            aniDict[99][nodeAni[4]] = [];
        }
        aniDict[99][nodeAni[4]].push(nodeHeight);*/
        }
        }
        // alternate way of getting branch heights matched to ANI ids, in case the above code for height fails
      /*$('text.tnt_tree_label:contains("ANIid_95")').each(function() {
        var aniIndex = $(this).text().search('ANIid_95');
        var aniId = $(this).text().slice(aniIndex+9).replace(/_/g, "");
        if (!(aniId in aniDict)) {
            aniDict[aniId] = [];
        }
        aniDict[aniId].push($(this).offset().top);
      });*/
      /*var aniSorted95 = sortAniList(aniDict[95]);
      var aniSorted96 = sortAniList(aniDict[96]);
      var aniSorted97 = sortAniList(aniDict[97]);
      var aniSorted98 = sortAniList(aniDict[98]);
      var aniSorted99 = sortAniList(aniDict[99]);*/
      var aniSorted = [sortAniList(aniDict[0]),sortAniList(aniDict[1]),sortAniList(aniDict[2]),sortAniList(aniDict[3]),sortAniList(aniDict[4])];
      //var prettyColors = ["#8B0000","#cc3600","#b36200","#7b7737","#556B2F","#006400"," #2d8653","#008080","#396a93","#00008B","#4B0082","#800080","#C71585","#DC143C"]; //sorted array of colors -> rainbow
      var prettyColors = [[ "#00990d", "#9900ff", " #e6007a", "#0000CD"],["#cc4400", "#008080", "#590099", "#518000"],["#e6007a", "#008080", "#590099", "#00990d"],["#cc4400", "#0000CD", "#9900ff", "#518000"],["#e6007a", "#0000CD", "#518000", "#590099"]];
      aniTable = [];
      aniTable[0] = {};
      aniTable[1] = {};
      aniTable[2] = {};
      aniTable[3] = {};
      aniTable[4] = {};
      // iterate through each list of ANIs, match up with corresponding color schemes
      for (var y in aniSorted) {
        for (var a in aniSorted[y]) {
        aniTable[y][aniSorted[y][a][1]] = prettyColors[y][a % prettyColors[y].length];
        }
      }
      // long version:
      /*aniTable[95] = {};
      aniTable[96] = {};
      aniTable[97] = {};
      aniTable[98] = {};
      aniTable[99] = {};*/
      /*for (var y in aniSorted95) {
      aniTable[95][aniSorted95[y][1]] = prettyColors[y % amountColors];
      }
      for (var a in aniSorted96) {
      aniTable[96][aniSorted96[a][1]] = ani96Colors[a % amountColors];
      }
      for (var b in aniSorted97) {
      aniTable[97][aniSorted97[b][1]] = ani97Colors[b % amountColors];
      }
      for (var c in aniSorted98) {
      aniTable[98][aniSorted98[c][1]] = ani98Colors[c % amountColors];
      }
      for (var d in aniSorted99) {
      aniTable[99][aniSorted99[d][1]] = ani99Colors[d % amountColors];
      }*/
      // coloring of single groups on click; "remembers" last ANI cutoff level via hidden input
        tree.on("click", function(node) {
            if (node.is_leaf()) {
            //var ani95 = node.property('ANI95');
            var lastAni = $('#lastani').val(); // values from ANI95-ANI99
            var currentNodeGroup = node.property(lastAni);
            //console.log(ani95);
            tree.label().color(function(node) {
                return subtreeFill(node,currentNodeGroup,lastAni);
            });
                        $('#colorkey').removeClass("hidden");
            $('#tscolorscheme').addClass('hidden');
            $('#anicolorscheme').addClass('hidden');
            $('#aniicon').css('color',subtreeFill(node,currentNodeGroup,lastAni));
            $('#selectedani').text(currentNodeGroup+'; cutoff: '+lastAni.slice(3)+'%');
            $('#singleani').removeClass('hidden');
            tree.update_nodes();
            } else {
                tree.label().color('black')
            }
            setSvg();
        });
        setSvg();
        // for printing; everything without .printable is hidden
        $('#svgCanvas').addClass('printable');
        $('#svgCanvas *').addClass('printable');

} else {
    $('#treefail').removeClass("hidden");
    $('#svgCanvas').addClass("hidden");
//    clearInterval(timer);
}
}

function treeError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}
function drawTree() {
$.ajax({
        url: '/results/'+jobid+'/tree',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: treeSuccess,
        error: treeError});

}

// tree resizing; only horizontal resizing possible
function resizeTree(resizeDir) {
    var treeWidth = tree.layout().width();

    if (resizeDir == 'upH') {
        tree.layout().width(treeWidth+100);
    } else if (resizeDir == 'downH' && treeWidth > 500) { //semi-arbitrary cutoff so display doesn't break
        tree.layout().width(treeWidth - 100);
    }
    tree.update();
    setSvg();
}
// returns font sizes used for search function
function treeFontSize(node,value) {
if (node.is_leaf()) {
    var name = node.node_name();
    if (name.toLowerCase().search(value.toLowerCase()) != -1) {
        return 14;
    }
    return 12;
}
    return 12;
    }
// searches tree; search results have larger font size
function searchTree() {
var value = $('#searchtree').val().toLowerCase();
if (value.length) {
    tree.label().fontsize(function(node) {
        return treeFontSize(node,value);
    });
    tree.update();
    }
    else {
    tree.label().fontsize(12);
    tree.update();
    }
}
// recolors tree according to ANI groups
function showAniGroups(aniCutoffLevel) {
if (aniCutoffLevel == '97') {
$('#colorblock1').css('color','#00990d');
$('#colorblock2').css('color','#9900ff');
$('#colorblock3').css('color','#e6007a');
$('#colorblock4').css('color','#0000CD');
$('#lastani').val('ANI97');
} else if (aniCutoffLevel == '95') {
$('#colorblock1').css('color','#cc4400');
$('#colorblock2').css('color','#008080');
$('#colorblock3').css('color','#590099');
$('#colorblock4').css('color','#518000');
$('#lastani').val('ANI95');
} else if (aniCutoffLevel == '90') {
$('#colorblock1').css('color','#e6007a');
$('#colorblock2').css('color','#008080');
$('#colorblock3').css('color','#590099');
$('#colorblock4').css('color','#00990d');
$('#lastani').val('ANI90');
} else if (aniCutoffLevel == '85') {
$('#colorblock1').css('color','#cc4400');
$('#colorblock2').css('color','#0000CD');
$('#colorblock3').css('color','#9900ff');
$('#colorblock4').css('color','#518000');
$('#lastani').val('ANI85');
} else if (aniCutoffLevel == '80') {
$('#colorblock1').css('color','#e6007a');
$('#colorblock2').css('color','#0000CD');
$('#colorblock3').css('color','#518000');
$('#colorblock4').css('color','#590099');
$('#lastani').val('ANI80');
}
$('#anicutoff').text(aniCutoffLevel);
$('#colorkey').removeClass('hidden');
$('#tscolorscheme').addClass('hidden');
$('#singleani').addClass('hidden');
$('#anicolorscheme').removeClass('hidden');
tree.label().color(function(node) {
if (node.is_leaf()) {
    return nodeFill(node,aniCutoffLevel);
    }
    else {
    return 'black';
    }
});
tree.update_nodes();
setSvg();
}
// recolors tree according to type strains/QS/outgroups
function showTypeStrains() {
$('#colorkey').removeClass('hidden');
$('#anicolorscheme').addClass('hidden');
$('#singleani').addClass('hidden');
$('#tscolorscheme').removeClass('hidden');
tree.label().color(function(node) {
    return strainTypeNodes(node);
});
tree.update_nodes();
setSvg();
}
// resets label colors to default (black)
function resetColor() {
$('#colorkey').addClass("hidden");
$('#tscolorscheme').addClass('hidden');
$('#anicolorscheme').addClass('hidden');
$('#singleani').addClass('hidden');
tree.label().color('black');
tree.update_nodes();
setSvg();
}

treeTimer = 0;

function startSearch() {
clearTimeout(treeTimer);
treeTimer = setTimeout(searchTree, 1000);
}

$(document).ready(function () {
$('#searchtree').on("keyup", startSearch);
});


// switches between scaled and unscaled layout
function toggleScale() {
if (tree.layout().scale()) {
    tree.layout().scale(false);
} else {
    tree.layout().scale(true);
    }
    tree.update();
    setSvg();
}
// opening tree SVG in new window if needed
function treePopup() {
    var treeWindow = window.open("","Tree View");
    treeWindow.document.write($('#svgCanvas *').html());
}

function aniError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

aniMasterList = null

function aniSuccess(data,textStatus,xhr) {
aniMasterList = data;
drawTree();
//console.log(aniMasterList);
}
// loading ANI groups for coloring
function loadAniGroups() {
$.ajax({
        // Your server script to process the upload
        url: '/static/aniclades.json',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: aniSuccess,
        error: aniError});
}


$("button[data-toggle='tooltip']").on('click',function(){
    $(this).blur();
})

$(document).ready(function() {
    loadAniGroups();
});