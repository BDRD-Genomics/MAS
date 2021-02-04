$(document).ready(function() {
    // COPIED FROM bdrd-lims.js
    t = $('.table').DataTable({
        pageLength: 10,
        lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
        // deferRender: true,
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

    t.column(':contains("Date Created")').order('desc').draw();
    t.columns.adjust();

    var feature = $('.featureTable').DataTable();
    var annotation = $('.annotationTable').DataTable();

    $('.selectmultiple:visible').select2({closeOnSelect: false});
    $('.select:visible').select2();

    feature.columns().every(function () {
        var that = this;

        $('input', this.footer()).on('keyup', function () {
            // console.log("feature table value: " + this.value);
            // console.log("feature table Search " + that.search());
            if (that.search() !== this.value) {
                that.search(this.value).draw();
            }
        });
    });

    annotation.columns().every(function () {
        var that = this;

        $('input', this.footer()).on('keyup', function () {
            // console.log("value: " + this.value);
            // console.log("Search " + that.search());
            if (that.search() !== this.value) {
                that.search(this.value).draw();
            }
        });
    });
});

$('div.genome-control').click(function() {
    var table_id = $(this).closest('table').attr('id');

    var row = t.table("#"+table_id).row($(this).closest('tr'));
    var phage_id = $(this).closest('tr').find("input:hidden[class='genome']").attr("value");

    if ( row.child.isShown() ) {
        // This row is already open - close it
        row.child.hide();
        $(this).removeClass('shown');
    }
    else{
        // Add class for CSS recognition
        $(this).addClass('shown');
        if (phage_id !== undefined) {
            $.ajax({
                url: "{% url 'genome:get_phage' %}",
                type: "GET",
                data: {'genome_id': phage_id},
                dataType: 'html',
                success: function (phage) {
                    row.child('<div id="phage_"'+phage_id+' class="phage">'+ phage +'</div>').show();
                }
            });
        }
    }

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

