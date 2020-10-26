function orgError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

function speciesSearch() {
    var searchval = $('#speciessearch').val().toLowerCase();
    $('#speciessel > option').each(function() {
        var test = ($(this).text().toLowerCase().indexOf(searchval) > -1);
        $(this).toggle(test);
        if (test && searchval.length) {
            $(this).css("color","#0000ff");
        } else {
            $(this).css("color","");
        }
    });
}

function outgrSearch() {
    var searchval = $('#outgroupsearch').val().toLowerCase();
    $('#outgrsel > option').each(function() {
        var test = ($(this).text().toLowerCase().indexOf(searchval) > -1);
        $(this).toggle(test);
        if (test && searchval.length) {
            $(this).css("color","#0000ff");
        } else {
            $(this).css("color","");
        }
        if (test && searchval.length && $(this).is(':disabled')) {
            $(this).toggle(false);
            $(this).css("color","");
        }
    });
}
var speciesTimer = 0;
//var speciesStarted = false;

/*function speciesQueue() {
    //console.log('SQ running');
    if (speciesTimer <= 0) {
        speciesSearch();
        speciesStarted = false;
    } else {
        speciesTimer -= 100;
        setTimeout(speciesQueue, 100);
    }
}*/
function speciesStart() {
clearTimeout(speciesTimer);
speciesTimer = setTimeout(speciesSearch, 1000);
}

$(document).ready(function() {
$('#speciessearch').on("keyup", speciesStart);
});

var outgrTimer = 0;
/*var outgrStarted = false;

function outgrQueue() {
    if (outgrTimer <= 0) {
        outgrSearch();
        outgrStarted = false;
    } else {
        outgrTimer -= 100;
        setTimeout(outgrQueue, 100);
    }
}*/

function outgrStart() {
clearTimeout(outgrTimer);
outgrTimer = setTimeout(outgrSearch, 1000);
}

$(document).ready(function() {
$('#outgroupsearch').on("keyup",outgrStart);
});
// add entries from object
function objectToSelectBox(inputObject, selectBox, isActive) {
    if (isActive == 'true') {
    $(selectBox).append("<option value='"+inputObject.id+"' data-title='"+ inputObject.orgname + " (mean distance: "+ parseFloat(inputObject.dist).toFixed(3)+")"+"' id='" + inputObject.id +"' data-genusid='"+inputObject.genusid+"' data-genusname='"+inputObject.genusname+"' data-familyid='"+inputObject.familyid+"' data-familyname='"+inputObject.familyname+"' data-orderid='"+inputObject.orderid + "' data-ordername='"+inputObject.ordername+"' data-phylid='"+inputObject.phylid+"' data-phylname='"+inputObject.phylname+"'>"+inputObject.orgname+ " (mean distance: "+ parseFloat(inputObject.dist).toFixed(3)+") </option>");
    } else {
    $(selectBox).append("<option value='"+inputObject.id+"' data-title='"+ inputObject.orgname + " (mean distance: "+ parseFloat(inputObject.dist).toFixed(3)+")"+"' id='" + inputObject.id +"' data-genusid='"+inputObject.genusid+"' data-genusname='"+inputObject.genusname+"' data-familyid='"+inputObject.familyid+"' data-familyname='"+inputObject.familyname+"' data-orderid='"+inputObject.orderid + "' data-ordername='"+inputObject.ordername+"' data-phylid='"+inputObject.phylid+"' data-phylname='"+inputObject.phylname+"' disabled>"+inputObject.orgname+ " (mean distance: "+ parseFloat(inputObject.dist).toFixed(3)+") </option>");
    }
}


var jsonOrgs = false;
// identify common taxonomic group of all selected organisms, as identified by NCBI ids
function commonGroup() {
var counter, counter2, commonTaxId;
//var selectedList = [];
var allGroup = {};
var groupList = ["genus","family","order","phyl"];
//var refOrgs = $.extend(jsonOrgs["reforgs"], jsonOrgs["queryorgs"]); //how to handle if only query orgs remaining?
 // for each taxonomical level, see what group each active option in species selection belongs to
 for (var group in groupList) {
                var groupId = groupList[group]+"id";
                var groupName = groupList[group]+"name";
                // iterate through all selected options
                $("#specieslist > .picked").each(function(){
                    var currId = $(this).val();
                    var currGroup = $(this).data(groupId);
                    var currName = $(this).data(groupName)
                    if (currGroup != "N/A" && currName != "N/A") {
                                allGroup[currGroup] = currName; //stores current option's group ID and name in allGroup; overwrites respective entry if already present
                                //console.log(allGroup);
                                }
                    /*for (counter=0; counter<refOrgs.length; counter++) {
                        var refInfo = refOrgs[counter];
                            if (currId == refInfo.id) {
                                //selectedList.push(refInfo);
                                var currGroup = refInfo[groupId];
                                var currName = refInfo[groupName];
                                if (currGroup != "N/A" && currName != "N/A") {
                                allGroup[currGroup] = currName;
                                console.log(allGroup);
                                }
                                }
                        }*/
                });
                // if all options belong to the same group at the level currently examined, return examined level, NCBI group ID and name of group
                if (Object.keys(allGroup).length == 1) {
                    commonTaxId = Object.keys(allGroup);
                    var commonInfo = [groupList[group], commonTaxId[0], allGroup[commonTaxId[0]]];
                //console.log(commonInfo);
                    return commonInfo;
                } else if (groupList[group] == "phyl" && Object.keys(allGroup).length > 1){ // if the check has been run down to phylum level without finding one single group to which all options belong, return level (phylum), as well as NCBI ids and names of all phyla represented
                    commonTaxId = Object.keys(allGroup);
                    var multiGroups = [];
                    for (var idkey in commonTaxId) {
                        multiGroups.push(allGroup[commonTaxId[idkey]]);
                    }
                    var commonInfo = [groupList[group], commonTaxId, multiGroups];
                    //console.log(commonInfo);
                    return commonInfo
                }
                // if neither of these conditions has been satisfied, empty allGroup and repeat check at next higher taxonomical level
                allGroup = {};
}
}

// select all organisms suitable as outgroups for current species selection, sort by Mash distance and add to outgroup list
function outgroupPick() {
var commonTax = commonGroup();
//console.log(typeof commonTax[1]);
var activeOutgroups = [];
var inactiveOutgroups = [];
if (typeof commonTax[1] == "string") { //if only one NCBI id in list
//var allOutgroups = jsonOrgs["outgroups"]; // change?
var groupId = commonTax[0]+"id";
$('#outgrsel').empty();
for (var counter in outgrInfo) { //iterate through outgroup dictionary
    if (outgrInfo[counter][groupId] == commonTax[1]) {
        inactiveOutgroups.push(outgrInfo[counter]);
        //console.log(outgrInfo[counter].id);
        removeFromList('#outgrlist',outgrInfo[counter].id);
    } else {
        activeOutgroups.push(outgrInfo[counter]);
    }

}
// sort outgroups from smallest to largest Mash distance
activeOutgroups.sort(function(a,b) {
    return a['dist'] - b['dist'];
});
inactiveOutgroups.sort(function(c,d) {
    return c['dist'] - d['dist'];
});
for (var counter3 in activeOutgroups) {
objectToSelectBox(activeOutgroups[counter3],'#outgrsel','true');
}
for (var counter4 in inactiveOutgroups) {
objectToSelectBox(inactiveOutgroups[counter4],'#outgrsel','false');
}
//var counter;
/*$('#outgrsel > option').each(function() {
    var currId = $(this).attr("id");
    if ($(this).data(groupId) == commonTax[1]) {
    $(this).attr("disabled", true);
                this.selected = false;
                removeFromList('#outgrlist', $(this).val());
    } else {
    $(this).removeAttr("disabled"); // faster to remove disabled from the start?
    }*/
    /*for (counter=0; counter<allOutgroups.length; counter++) {
        var outgrInfo = allOutgroups[counter];
        if (currId == outgrInfo.id) {
            if (outgrInfo[groupId] == commonTax[1]) {
                console.log(outgrInfo[groupId]);
                $(this).attr("disabled", true);
                this.selected = false;
                removeFromList('#outgrlist', $(this).val());
            } else {
                $(this).removeAttr("disabled");
            }
        }
    }
});*/
$('#commgroup').text(commonTax[2]);
if (!($('#outgrsel > option:not([disabled])').length)) {
    $('#nooutgroups').removeClass("hidden");
} else {
    $('#nooutgroups').addClass("hidden");
}
refreshView('#outgrlist',5,'#outgrlimit');
outgrSearch();
} else if (typeof commonTax[1] == "object") { //if more than one NCBI id in group
$('#outgrsel').empty();
    //var allOutgroups = jsonOrgs["outgroups"]; // change?
    var groupId = commonTax[0]+"id";
    //var counter2;
    for (var counter2 in outgrInfo) {
        if (outgrInfo[counter2][groupId] in commonTax[1]) {
            inactiveOutgroups.push(outgrInfo[counter2]);
            removeFromList('#outgrlist',outgrInfo[counter2].id);
        } else {
            inactiveOutgroups.push(outgrInfo[counter2]);
        }
    }
    // sort outgroups from smallest to largest Mash distance
    activeOutgroups.sort(function(e,f) {
        return e['dist'] - f['dist'];
    });
    inactiveOutgroups.sort(function(g,h) {
        return g['dist'] - h['dist'];
    });
    for (var counter5 in activeOutgroups) {
        objectToSelectBox(activeOutgroups[counter5],'#outgrsel','true');
    }
    for (var counter6 in inactiveOutgroups) {
        objectToSelectBox(inactiveOutgroups[counter6],'#outgrsel','false');
    }
    /*$('#outgrsel > option').each(function() {
    //var currId = $(this).attr("id");
    $(this).removeAttr("disabled");
    for (counter2=0; counter2<commonTax[1].length; counter2++) {
            currTax = commonTax[1][counter2];
            if ($(this).data(groupId) == currTax) {
            $(this).attr("disabled", true);
            this.selected = false;
            removeFromList('#outgrlist', $(this).val());
            $(this).remove();
            $('#outgrsel').append($(this));
            console.log(currTax);
            }
            }
    /*for (counter=0; counter<allOutgroups.length; counter++) {
    var outgrInfo = allOutgroups[counter];
    if (currId == outgrInfo.id) {
        for (counter2=0; counter2<commonTax[1].length; counter2++) {
            currTax = commonTax[1][counter2];
            if (outgrInfo[groupId] == currTax) {
            $(this).attr("disabled", true);
            this.selected = false;
            removeFromList('#outgrlist', $(this).val());
            $(this).remove();
            $('#outgrsel').append($(this));
            console.log(currTax);
            }
        }
    }
    }
    });*/
    $('#commgroup').text(commonTax[2].join(", "));
if (!($('#outgrsel > option:not([disabled])').length)) {
    $('#nooutgroups').removeClass("hidden");
} else {
    $('#nooutgroups').addClass("hidden");
}
refreshView('#outgrlist',5,'#outgrlimit');
outgrSearch();
} else {
    $('#nooutgroups').removeClass("hidden");
    $('#outgrsel > option').each(function() {
    $(this).attr("disabled", true);
    this.selected = false;
    removeFromList('#outgrlist', $(this).val());
    $('#commgroup').text("No common group found");
    });
    outgrSearch();
}
}

// add all query and reference organisms to species multiselect; set outgroup information as global variable, fill in from data, then add to outgroup multiselect
function orgSuccess(data,textStatus,xhr) {
var counter;
var counter2;
var counter3;
jsonOrgs = data;
var queryOrgs = data["queryorgs"];
var refOrgs = data["reforgs"];
var outgroups = data["outgroups"];
outgrInfo = [];
for (counter3 = 0; counter3<queryOrgs.length; counter3++) {
    var queryInfo = queryOrgs[counter3];
    $('#speciessel').append("<option value='"+queryInfo.id+"' data-title='(U) "+ queryInfo.orgname + " (detected genus: "+ queryInfo.genusname +")"+"' id='" + queryInfo.id +"' data-genusid='"+queryInfo.genusid+"' data-genusname='"+queryInfo.genusname+"' data-familyid='"+queryInfo.familyid+"' data-familyname='"+queryInfo.familyname+"' data-orderid='"+queryInfo.orderid + "' data-ordername='"+queryInfo.ordername+"' data-phylid='"+queryInfo.phylid+"' data-phylname='"+queryInfo.phylname+"'> (U) "+queryInfo.orgname+ " (detected genus: "+ queryInfo.genusname+") </option>");
}
for (counter = 0; counter<refOrgs.length; counter++) {
    var orgInfo = refOrgs[counter];
    typestr = (orgInfo.typestrain) ? orgInfo.typestrain : "false";
    if (typestr == true) {
    $('#speciessel').append("<option value='"+orgInfo.id+"' data-title='(T) "+ orgInfo.orgname + " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+")"+"' id='" + orgInfo.id + "' data-genusid='"+orgInfo.genusid+"' data-genusname='"+orgInfo.genusname+"' data-familyid='"+orgInfo.familyid+"' data-familyname='"+orgInfo.familyname+"' data-orderid='"+orgInfo.orderid + "' data-ordername='"+orgInfo.ordername+"' data-phylid='"+orgInfo.phylid+"' data-phylname='"+orgInfo.phylname+"'>(T) "+orgInfo.orgname+ " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+") </option>");
    } else {
    $('#speciessel').append("<option value='"+orgInfo.id+"' data-title='"+ orgInfo.orgname + " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+")"+"' id='" + orgInfo.id +"' data-genusid='"+orgInfo.genusid+"' data-genusname='"+orgInfo.genusname+"' data-familyid='"+orgInfo.familyid+"' data-familyname='"+orgInfo.familyname+"' data-orderid='"+orgInfo.orderid + "' data-ordername='"+orgInfo.ordername+"' data-phylid='"+orgInfo.phylid+"' data-phylname='"+orgInfo.phylname+"'>"+orgInfo.orgname+ " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+") </option>");
   }
   }
for (counter2 = 0; counter2<outgroups.length; counter2++) {
    outgrInfo[counter2] = outgroups[counter2]; //add entries to object; outgroup multiselect is later populated by outgroupPick
    //$('#outgrsel').append("<option value='"+outgrInfo[counter2].id+"' data-title='"+ outgrInfo[counter2].orgname + " (mean distance: "+ parseFloat(outgrInfo[counter2].dist).toFixed(3)+")"+"' id='" + outgrInfo[counter2].id +"' data-genusid='"+outgrInfo.genusid+"' data-genusname='"+outgrInfo.genusname+"' data-familyid='"+outgrInfo.familyid+"' data-familyname='"+outgrInfo[counter2].familyname+"' data-orderid='"+outgrInfo[counter2].orderid + "' data-ordername='"+outgrInfo[counter2].ordername+"' data-phylid='"+outgrInfo[counter2].phylid+"' data-phylname='"+outgrInfo[counter2].phylname+"'>"+outgrInfo[counter2].orgname+ " (mean distance: "+ parseFloat(outgrInfo[counter2].dist).toFixed(3)+") </option>");

}
if (data["selspecies"] && data["seloutgroups"]) {
    var selectedSpecies = data["selspecies"];
    var selectedOutgroups = data["seloutgroups"];
    loadSelected('#speciessel','#specieslist',selectedSpecies);
    outgroupPick();
    loadSelected('#outgrsel','#outgrlist',selectedOutgroups);
} else {
    loadDefaults('#speciessel','#specieslist',50);
    outgroupPick();
    loadDefaults('#outgrsel','#outgrlist',1);
}
}
// load initial species and outgroups
function getOrgs() {
var jobid = $('#jobinfo').val();
//console.log(jobid)
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results/'+jobid+'/step2/orgs',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: orgSuccess,
        error: orgError});

}

function outgroupError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

function outgroupSuccess(data,textStatus,xhr) {
//console.log(data);
var counter;
//var counter2;
//    var optionList = [];
for (counter=0; counter<data.length;counter++) { //if setup with outgrInfo works, only needs to add to outgrInfo, the rest is handled by outgroupPick
    //var newOutgrInfo[counter] = data[counter];
    if (!(data[counter] in outgrInfo)) {
    outgrInfo.push(data[counter]);
    }
    /*$('#outgrsel > option').remove('#'+newOutgrInfo.id); //fix ordering?
    var optionInfo = "<option value='"+newOutgrInfo.id+"' data-title='"+ newOutgrInfo.orgname + " (mean distance: "+ parseFloat(newOutgrInfo.dist).toFixed(3)+")"+"' id='" + newOutgrInfo.id +"' data-genusid='"+newOutgrInfo.genusid+"' data-genusname='"+outgrInfo.genusname+"' data-familyid='"+newOutgrInfo.familyid+"' data-familyname='"+newOutgrInfo.familyname+"' data-orderid='"+newOutgrInfo.orderid + "' data-ordername='"+newOutgrInfo.ordername+"' data-phylid='"+newOutgrInfo.phylid+"' data-phylname='"+newOutgrInfo.phylname+"'>"+newOutgrInfo.orgname+ " (mean distance: "+ parseFloat(newOutgrInfo.dist).toFixed(3)+") </option>";
    optionList.unshift(optionInfo);
}
for (counter2=0; counter2<optionList.length; counter2++) {
    $('#outgrsel').prepend(optionList[counter2]);*/
}
outgroupPick(); //handles sorting by Mash distance, filtering of outgroups, adding outgroups to multiselect
if (!($('#outgrsel > option:not([disabled])').length)) {
    $('#nooutgroups').removeClass("hidden");
} else {
    $('#nooutgroups').addClass("hidden");
}
$('#outgroupsearch').val("");
outgrSearch();
}
// loading more outgroups, pre-filtered to exclude any belonging to the selected organisms' common  taxonomical group
function outgroupLoad() {
var jobid = $('#jobinfo').val();
var commonTax = commonGroup();
//add check for if no options in outgrsel?
if (typeof commonTax[1] == "string") {
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results/'+jobid+'/step2/outgroups?group='+commonTax[0]+'&id='+commonTax[1],
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: outgroupSuccess,
        error: outgroupError});
} else if (typeof commonTax[1] == "object") {
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results/'+jobid+'/step2/outgroups?group='+commonTax[0]+'&id='+commonTax[1].toString()+'&multiple=True',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: outgroupSuccess,
        error: outgroupError});
}
}


function moreSeqsError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

function moreSeqsSuccess(data,textStatus,xhr) {
var counter;
var refOrgs=data["reforgs"];
jsonOrgs["reforgs"] = $.extend(jsonOrgs["reforgs"],data["reforgs"]);
for (counter=0; counter<refOrgs.length;counter++) {
    var orgInfo = refOrgs[counter];
    typestr = (orgInfo.typestrain) ? orgInfo.typestrain : "false";
    if (typestr == true) {
    $('#speciessel').append("<option value='"+orgInfo.id+"' data-title='(T) "+ orgInfo.orgname + " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+")"+"' id='" + orgInfo.id +"' data-genusid='"+orgInfo.genusid+"' data-genusname='"+orgInfo.genusname+"' data-familyid='"+orgInfo.familyid+"' data-familyname='"+orgInfo.familyname+"' data-orderid='"+orgInfo.orderid + "' data-ordername='"+orgInfo.ordername+"' data-phylid='"+orgInfo.phylid+"' data-phylname='"+orgInfo.phylname+"'>(T) "+orgInfo.orgname+ " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+") </option>");
    } else {
    $('#speciessel').append("<option value='"+orgInfo.id+"' data-title='"+ orgInfo.orgname + " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+")"+"' id='" + orgInfo.id +"' data-genusid='"+orgInfo.genusid+"' data-genusname='"+orgInfo.genusname+"' data-familyid='"+orgInfo.familyid+"' data-familyname='"+orgInfo.familyname+"' data-orderid='"+orgInfo.orderid + "' data-ordername='"+orgInfo.ordername+"' data-phylid='"+orgInfo.phylid+"' data-phylname='"+orgInfo.phylname+"'>"+orgInfo.orgname+ " (mean distance: "+ parseFloat(orgInfo.dist).toFixed(3)+") </option>");
   }
    }
    refInfo.rows.add(data["reforgs"]).draw();
    $('#speciessearch').val("");
    speciesSearch();
}

// loading more sequences for species multiselect
function loadMore() {
var jobid = $('#jobinfo').val();
var loadedSeqs = $('#speciessel > option').length - jsonOrgs["queryorgs"].length; // amount of species loaded, minus user-submitted species
//console.log(loadedSeqs);
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results/'+jobid+'/step2/orgs?start='+loadedSeqs,
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: moreSeqsSuccess,
        error: moreSeqsError});
}



/*$(document).ready(function() {
    $('#speciessearch').on("keyup", function() {
        var searchval = $(this).val().toLowerCase();
        $('#speciessel > option').each(function() {
        var test = ($(this).text().toLowerCase().indexOf(searchval) > -1);
        $(this).toggle(test);
        if (test && searchval.length) {
            $(this).css("color","#0000ff");
        } else {
            $(this).css("color","");
        }
        });
    });
});

$(document).ready(function() {
    $('#outgroupsearch').on("keyup", function() {
        var searchval2 = $(this).val().toLowerCase();
        $('#outgrsel > option').each(function() {
        var test = ($(this).text().toLowerCase().indexOf(searchval2) > -1);
        $(this).toggle(test);
        if (test && searchval2.length) {
            $(this).css("color","#0000ff");
        } else {
            $(this).css("color","");
        }
        });
    });
}); */

// check if enough organisms (at least 5) and enough outgroups (at least 1) selected
function validateForm() {
if ($('.selectablein >option').length>4 && $('.selectablein2 > option').length>0) {
return "validated";
} else {
return "notvalidated";}
}