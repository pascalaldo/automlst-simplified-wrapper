var genusKeys = false;
// lists genera in to which user sequences belong; if there are multiple, user needs to select one to proceed with
// may need some troubleshooting
function genusSuccess(data,textStatus,xhr) {
    var genusList = data["genuslist"];
    if (String(genusKeys) != String(Object.keys(genusList))) {
    genusKeys = Object.keys(genusList);
    //need to clear the counter?
    //console.log(genusKeys.length);
    /*var genusVals = Object.values(genusList);
    var maxValue = Math.max.apply(null,genusVals); */
    var queries = data["queryfiles"];
    var maxGenus = data["maxgenus"];
    //console.log(maxGenus);
    var counter;
    var counter2;
    // if more than one genus represented, list of genera to select from built; table with information built
    if (genusKeys.length > 1) {
        $('#genuswarn').removeClass("hidden");
        $('#genustable >tbody').empty();
        for (counter=0;counter<genusKeys.length;counter++) {
            $('#genusoptions').append("<option id='"+genusKeys[counter]+"' value='"+genusKeys[counter]+"'>"+genusKeys[counter]+": "+genusList[genusKeys[counter]]+" sequences </option>");
            if (genusKeys[counter] == maxGenus) {
                //console.log(genusKeys[counter]);
                $('#'+genusKeys[counter]).attr("selected","selected");
            }
        }
        for (counter2=0; counter2<queries.length;counter2++) {
            var currentSeq = queries[counter2];
            $('#genustable >tbody').append("<tr><td>"+currentSeq.file+"</td><td>"+currentSeq.genus+"</td></tr>");
        }
    } else if (genusKeys.length == 1) {
        $('#selectedgenus').text(maxGenus);
        $('#genuswarn').addClass("hidden");
    } else {
        $('#genuswarning').removeClass("hidden");
    }
}
}

function genusError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
}

function getRefGenus() {
var jobid = $('#jobinfo').val();
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results2/'+jobid+'/refs',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: genusSuccess,
        error: genusError});
}
function selectError(xhr,ajaxOptions,thrownError) { // communication error
console.log(xhr,ajaxOptions,thrownError);
$('#communicateerror').removeClass("hidden");
}
function selectSuccess(data,textStatus,xhr){
var data = JSON.parse(data);
if (parseInt(data["status"]) != 1) {
    $('#statuserror').removeClass("hidden");
} else {
//$('#selectedgenus').text(genusKeys[0]);
$('#selectedgenus').text(data["newmax"]);
        $('#genuswarn').addClass("hidden");
}
}
function selectGenus() {
$('#communicateerror').addClass("hidden");
$.ajax({
        // Your server script to process the upload
        contentType: 'application/json',
        url: '/results2/selectgenus',
        async: true,
        type: 'POST',
        //Form data
        data: new FormData($("#multigenus")[0]),
        cache: false,
        contentType: false,
        processData: false,
        success: selectSuccess,
        error: selectError});
}