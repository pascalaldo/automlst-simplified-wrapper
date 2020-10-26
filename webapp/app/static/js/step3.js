// displays organisms that need to be deleted to proceed with current gene selection
function findDeleted() {
var deletedList = [];
$('#deletewarning').addClass('hidden');
$('#todelete').empty();
$('#mlstlist > option').each(function() { // which organisms don't contain the gene in question?
    if ($(this).data("del")) { //if gene requires organisms to be deleted (recorded in data-del)
        $(this).css({"background-color":"#fcf8e5"}); // highlights background if organisms need to be deleted
        var currDel = $(this).data("del").split(",");
        var i;
        for (i in currDel) {
            if (deletedList.indexOf(currDel[i]) == -1) { // add to list of deleted organisms if not already in the list
                deletedList.push(currDel[i]);
            }
        }
    }
});
if (deletedList.length > 0) { // if there are deleted organisms, show list
    var j;
    for (j in deletedList) {
        $('#todelete').append("<li>"+deletedList[j]+"</li>"); // add entries to list displayed to user
    }
    $('#deletewarning').removeClass("hidden"); // *then* show the user that list
    $('#removeorgs').val(deletedList.toString()); // currently without function...
}
}

function geneError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

// load genes into multiselects, add information
function geneSuccess(data,textStatus,xhr) {
    var counter;
    var mlstGenes = data;
    /*var mlstGenes = data["genes"];
    var selectedGenes = data["selected"];*/
    for (counter=0; counter<mlstGenes.length; counter++) {
        var geneInfo = mlstGenes[counter];
        //console.log(geneInfo);
        if (geneInfo["orgdel"].length>0) {
        $('#mlstsel').append("<option id='" + geneInfo.acc + "' data-title='"+ geneInfo.name +": "+geneInfo.desc+"' value='"+geneInfo.acc+"' data-del='"+geneInfo.orgdel +"' style='background-color:#fcf8e5'>" + geneInfo.name +": "+geneInfo.desc+"</option>");
        } else {
        $('#mlstsel').append("<option id='" + geneInfo.acc + "' data-title='"+ geneInfo.name +": "+geneInfo.desc+"' value='"+geneInfo.acc+"' data-del='"+geneInfo.orgdel +"'>" + geneInfo.name +": "+geneInfo.desc+"</option>");
        }
        }
     if (data["selected"]) { //currently an useless check
        loadSelected('#mlstsel','#mlstlist',selectedGenes);
     } else {
        loadDefaults('#mlstsel','#mlstlist',100);
    }
    findDeleted();
    /*if (data["deletedorgs"]) {
        var counter2;
        for (counter2=0;counter2<deletedOrgs.length;counter2++){
            var delInfo = deletedOrgs[counter2];
            $('#todelete').append("<li>"+delInfo.orgname+"</li>");
        }
        $('#deletewarning').removeClass("hidden");
    } */
}

function getGenes() {
var jobid = $('#jobinfo').val();
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results/'+jobid+'/step3/genes',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: geneSuccess,
        error: geneError});

}


function validateForm() {
if ($('.selectablein>option').length>10) {
return "validated";
} else {
return "notvalidated";}
}