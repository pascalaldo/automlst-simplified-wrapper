var jobid = $('#jobinfo').val();
$(document).ready(function() {
     mashTable = $('#mashinfo').DataTable({
        ajax: {
         contentType: 'application/json',
         url:'/results/'+jobid+'/mash'
        },
        columns: [
            {data : 0,
            responsivePriority:1},
            {data: 1,
            responsivePriority:6},
            {data: 2,
            responsivePriority: 2,
            render: function(data,type,row,meta){
                return data.replace(/_/g," ");
            }},
            {data: 3,
            responsivePriority: 5},
            {data: 4,
            responsivePriority: 3,
            render: function(data,type,row,meta){
                return (data*100).toFixed(1) + "%";}},
            {data: 5,
            responsivePriority: 9},
            {data: 6,
            responsivePriority: 7},
            {data: 7,
            responsivePriority: 8},
            {data: 8,
            responsivePriority: 4}
        ],
        order: [[4, 'desc']],
    responsive: {
        details: {
        display: $.fn.dataTable.Responsive.display.childRow
                }
    }

    });

    reloadTimer = setInterval( function () {
        if ( ! mashTable.data().count() ){
            mashTable.ajax.reload( null, false );
        } else {
            clearInterval(reloadTimer);
        }
    }, 3000 );
});