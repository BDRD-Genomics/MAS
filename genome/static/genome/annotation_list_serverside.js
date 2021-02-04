$(document).ready(function() {
    var AnnotationTable = $('#annotation-table').DataTable({
        order: [[ 2, "asc" ]],
        orderClasses: false,
        serverSide: true,
        processing: true,
        ajax: {
            url: GET_ANNO_DATA_URL,
            type: "GET",
            "data": {genome_id: GENOME_ID}
        },
        columns: [
            { "data": "sequence" },
            { "data": "feature" },
            { "data": "accession" },
            { "data": "annotation" },
            { "data": "length" },
            { "data": "public_notes" },
            { "data": "private_notes" },
            { "data": "flag" },
            { "data": "assigned_to" },
            { "data": "download" },
            { "data": "history" },
            { "data": "view_results" }
        ],
        columnDefs: [
            {
                targets: [0, 1, 9, 10, 11], /* column index */
                orderable: false, /* true or false */
            },
            { className: "centered", targets: [ 9, 10, 11 ] }
        ]
    });
    $('table').on('click', 'tbody > tr > td > div.aa-control', function() {
        var table_id = $(this).closest('table').attr('id');

        var row = AnnotationTable.table("#"+table_id).row($(this).closest('tr'));

        var annotation_id = $(this).closest('td').find("input:hidden[class='annotation']").attr("value");

        if ( row.child.isShown() )
        {
            // This row is already open - close it
            row.child.hide();
            $(this).removeClass('shown');
        }
        else{
            // Add class for CSS recognition
            $(this).addClass('shown');
            if (annotation_id !== undefined) {
                $.ajax({
                    // url: "{% url 'genome:get_aa_sequence' %}",
                    url: GET_AA_SEQ_URL,
                    type: "GET",
                    data: {'annotation_id': annotation_id},
                    dataType: 'html',
                    success: function (annotation) {
                        row.child('<div id="annotation_"'+annotation_id+' class="annotation">'+ annotation +'</div>').show();
                    }
                });
            }
        }
    });
});