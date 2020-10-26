var jobid = $('#jobinfo').val();
$(document).ready(function() {
mlstInfo = $('#geneinfo').DataTable({
    ajax: {
        contentType: 'application/json',
        url: '/results/'+jobid+'/step3/genes',
        //dataSrc: 'genes'
        dataSrc: ''
    },
    columns: [
        {data: "name",
        responsivePriority: 1},
        {data: "acc",
        responsivePriority: 2},
        {data: "func",
        responsivePriority: 4},
        {data: "desc",
        responsivePriority: 3}
    ],
    order: [[0, "asc"]],
        responsive: {
        details: {
        display: $.fn.dataTable.Responsive.display.childRow
                }
    },
    aLengthMenu: [
        [10, 25, 50, 100, 200, -1],
        [10, 25, 50, 100, 200, "All"]
        ],
        iDisplayLength: 10,
     fnInitComplete: function(oSettings, json) {
      idSearch();
    }
    });
 });

function idSearch() {
var idList = "";
$('#mlstlist .picked').each(function(){
    idList += ("|"+($(this).val()));
});
//console.log(idList);
var searchTerm = idList.slice(1);
mlstInfo.search(searchTerm, true,false).draw();
}