//tried to add column search
//need documentation to do more
$(document).ready(function() {



    var entireTable = $('#table').dataTable();

    var allpages = entireTable.fnGetNodes();

    $('#select-all').click(function () {
        if ($(this).hasClass('allChecked')){
            $('input[type="checkbox"]', allpages).prop('checked',false);
             $(this).toggleClass('allChecked')
        }
        else{
            $('input[type="checkbox"]', allpages).prop('checked',true);
            $(this).toggleClass('allChecked')
        }
    });
});

// $('div.aa-control').click(function() {
//      console.log('AA - I have been clicked');
//     var table_id = $(this).closest('table').attr('id');
//
//     var row = t.table("#"+table_id).row($(this).closest('tr'));
//     // console.log("table id: " + table_id);
//     // console.log("row:" + row);
//     var annotation_id = $(this).closest('tr').find("input:hidden[class='annotation']").attr("value");
//     // console.log("phage: " + annotation_id);
//
//     if ( row.child.isShown() ) {
//         // This row is already open - close it
//         row.child.hide();
//         $(this).removeClass('shown');
//     }
//     else{
//         // Add class for CSS recognition
//         $(this).addClass('shown');
//         if (annotation_id !== undefined) {
//             $.ajax({
//                 // url: "{% url 'sop:get_qc_info' %}",
//                 url: "{% url 'genome:get_aa_sequence' %}",
//                 // url: "/genome/ajax/get_aa_sequence",
//                 type: "GET",
//                 data: {'annotation_id': annotation_id},
//                 dataType: 'html',
//                 success: function (annotation) {
//                     row.child('<div id="annotation_"'+annotation_id+' class="annotation">'+ annotation +'</div>').show();
//                 }
//             });
//         }
//     }
// });




