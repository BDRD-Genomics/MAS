$(document).ready(function() {
    var feature_table = $('.featureTable').DataTable({
        pageLength: 10,
        lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
        scrollY: "100%",
        scrollX: true,
        footer: true,
        buttons: [
            {
                extend: 'copyHtml5',
                text: 'Copy All'
            },
            {
                extend: 'copyHtml5',
                text: 'Copy Visible',
                exportOptions: {columns: ':visible'}
            },
            {
                extend: 'excelHtml5',
                text: 'Excel All'
            },
            {
                extend: 'excelHtml5',
                text: 'Excel Visible',
                exportOptions: {columns: ':visible'}
            },
            {
                extend: 'csvHtml5',
                text: 'TSV All',
                fieldSeparator: '\t', extension: '.txt',
                fieldBoundary: ''
            },
            {
                extend: 'csvHtml5',
                text: 'TSV Visible',
                fieldSeparator: '\t',
                extension: '.txt',
                fieldBoundary: '',
                exportOptions: {columns: ':visible'}
            }
        ],
    });

    $('table').on('click', 'tbody > tr > td > div.sequence-control', function () {
        console.log('seqcon - I have been clicked');
        var table_id = $(this).closest('table').attr('id');

        var row = feature_table.table("#"+table_id).row($(this).closest('tr'));
        var feature_id = $(this).closest('tr').find("input:hidden[class='feature']").attr("value");

        if ( row.child.isShown() ) {
            // This row is already open - close it
            row.child.hide();
            $(this).removeClass('shown');
        }
        else{
            // Add class for CSS recognition
            $(this).addClass('shown');
            if (feature_id !== undefined) {
                $.ajax({
                    url: GET_FEATURE_SEQ_URL,
                    type: "GET",
                    data: {'feature_id': feature_id},
                    dataType: 'html',
                    success: function (feature) {
                        row.child('<div id="feature_"'+feature_id+' class="feature">'+ feature +'</div>').show();
                    }
                });
            }
        }
    });
});