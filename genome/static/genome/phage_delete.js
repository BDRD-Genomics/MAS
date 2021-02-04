// $('button.button').click(function() {
//     var phage_id_list = [];
//     // var table = documnet.getElementById("table");
//     // for(var i = 0, row; row = table.rows[i]; i++){
//     //     for(var j = 0, col; col = row.cells[j]; j++){
//     //
//     //     }
//     // }
//     $('#table').find('input[type="checkbox"]:checked').each(function () {
//         var phage_id = $(this).closest('tr').find("input:checkbox[class='action-select']").attr("value");
//         console.log("phage id: " + phage_id);
//         phage_id_list.push(phage_id);
//     });
//
//     console.log("phage id list: " + phage_id_list);
//
//     $.ajax({
//         // url: "{% url 'sop:get_qc_info' %}",
//         url: "/genome/ajax/phage/delete",
//         type: "GET",
//         data: {'phage_id_list': phage_id_list},
//         dataType: 'html',
//         // return information of template
//         // success: function (phage_id) {
//         //     // row.child('<div id="annotation_"'+annotation_id+' class="annotation">'+ annotation +'</div>').show();
//         //     ('<div> success </div>').show();
//         // }
//
//     });
//
// });