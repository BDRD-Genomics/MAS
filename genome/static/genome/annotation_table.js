
$(document).ready(function() {
    // From https://www.gyrocode.com/articles/jquery-datatables-how-to-submit-all-pages-form-data/
    $("form").on('submit', function(e){

        var form = this;
        var params = t.$('input, select').serializeArray();

        $.each(params, function(){
            // If element doesn't exist in DOM
            if(!$.contains(document, form[this.name])){
            // Create a hidden element
                $(form).append(
                    $('<input>')
                      .attr('type', 'hidden')
                      .attr('name', this.name)
                      .val(this.value)
                );
            }
        });
    } );

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
//      //console.log('I have been clicked');
//     console.log('showing')
//     var table_id = $(this).closest('table').attr('id');
//
//     var row = t.table("#"+table_id).row($(this).closest('tr'));
//
//     var annotation_id = $(this).closest('tr').find("input:hidden[class='annotation']").attr("value");
//
//     if ( row.child.isShown() )
//     {
//         // This row is already open - close it
//         row.child.hide();
//         $(this).removeClass('shown');
//     }
//     else{
//         // Add class for CSS recognition
//         $(this).addClass('shown');
//         if (annotation_id !== undefined) {
//             $.ajax({
//                 // url: "{% url 'genome:get_aa_sequence' %}",
//                 url: GET_AA_SEQ_URL,
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

// $('table').on('click', 'tbody > tr > td > div.aa-control', function() {
//         var table_id = $(this).closest('table').attr('id');
//
//         var row = t.table("#"+table_id).row($(this).closest('tr'));
//
//         var annotation_id = $(this).closest('td').find("input:hidden[class='annotation']").attr("value");
//
//         if ( row.child.isShown() )
//         {
//             // This row is already open - close it
//             row.child.hide();
//             $(this).removeClass('shown');
//         }
//         else{
//             // Add class for CSS recognition
//             $(this).addClass('shown');
//             if (annotation_id !== undefined) {
//                 $.ajax({
//                     // url: "{% url 'genome:get_aa_sequence' %}",
//                     url: GET_AA_SEQ_URL,
//                     type: "GET",
//                     data: {'annotation_id': annotation_id},
//                     dataType: 'html',
//                     success: function (annotation) {
//                         row.child('<div id="annotation_"'+annotation_id+' class="annotation">'+ annotation +'</div>').show();
//                     }
//                 });
//             }
//         }
//     });
