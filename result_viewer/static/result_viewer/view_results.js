
var BAR_HEIGHT = 18;
var CHART_WIDTH = 800;
var CONTAINER_WIDTH = 850;
var KEY_SPACING = 80;
var ROW_SPACING = 10;

var container = d3.select('div#viz')
    .style('width', `${CONTAINER_WIDTH}px`);

var s = d3.scaleLinear()
    .domain([0, QUERY_LEN])
    .range([0, CHART_WIDTH]);

var xAxis = d3.axisBottom()
    .scale(s);

var annotation_in_db = JSON.parse(document.getElementById("current_annotation_info").textContent);

// Create color scale for key and subjects
var colorScale;  // initialized in update()

var alignment = d3.select("#alignment")
    .attr("width", CHART_WIDTH)
    .attr("height",  10)
    .on('mousemove', function () {
        var tick_locs = xAxis.tickValues();
        tick_locs.pop();
        tick_locs.push(Math.round(s.invert(d3.mouse(this)[0])));
        xAxis.tickValues(tick_locs);
        d3.select('#query').call(xAxis);
        d3.selectAll('#query g.tick').filter(function (d, i, list) {
            return i === list.length - 1;
        }).classed('mouse_tick', true);
    });


function display_correct_db_buttons(button_group_id) {
    $('#database-button-div').children('.button-group').each(function () {
        var group = $(this);

        if (group.attr('id') === button_group_id) {
            group.show();
            group.children('button').each(function () {
                var button = $(this);

                if (button.text() === DATABASE) {
                    button.attr('disabled', 'disabled');

                } else {
                    button.attr('disabled', null);
                }
            });

        } else {
            group.hide();
        }
    });
}


function on_tool_button_click(tool) {
    // Switch tool
    TOOL = tool;

    // Enable other buttons - Disable this button
    $('#tool-selection-buttons').children('button').each(function () {
        var button = $(this);

        if (button.attr('id') === `${tool}_button`) {
            button.attr('disabled', 'disabled');

        } else {
            button.attr('disabled', null);
        }
    });

    if (TOOL === 'blastp') {
        DATABASE = Object.keys(blastp_alignment_data)[0];
    } else if (TOOL === 'hhsearch') {
        DATABASE = Object.keys(hhsearch_alignment_data)[0];
    } else if (TOOL === 'rpsblast') {
        DATABASE = Object.keys(rpsblast_alignment_data)[0];
    }

    display_correct_db_buttons(`${tool}-db-selection-buttons`);

    // Remove selected class from table
    $('tr.selected').removeClass('selected');

    update();
}


d3.select('#hhsearch_button')
    .on('click', function () {
        if (TOOL !== 'hhsearch') {
            on_tool_button_click('hhsearch');;
        }
    });


d3.select('#blastp_button')
    .on('click', function () {
        if (TOOL !== 'blastp') {
            on_tool_button_click('blastp')
        }
    });


d3.select('#rpsblast_button')
    .on('click', function () {
        if (TOOL !== 'rpsblast') {
            on_tool_button_click('rpsblast')
        }
    });


$('.db-button').each(function () {
    var button = $(this)

    d3.select(this).on('click', function () {
        DATABASE = button.text();

        $('.db-button').each(function () {
            $(this).attr('disabled', null);
        });
        button.attr('disabled', 'disabled');

        update();
    });
});


d3.select('#svg-row-selector')
    .on('change', function () {
        update();
    });


d3.select('#copy-sequence')
    .on('click', function () {
        console.log('Clicked***');
       var temp_el = document.createElement('textarea');
       // Set value (string to be copied)
       temp_el.value = SEQ;
       // Set non-editable to avoid focus and move outside of view
       temp_el.setAttribute('readonly', '');
       temp_el.style = {position: 'absolute', left: '-9999px'};
       document.body.appendChild(temp_el);
       // Select text inside element
       temp_el.select();
       // Copy text to clipboard
       document.execCommand('copy');
       // Remove temporary element
       document.body.removeChild(temp_el);
    });


d3.select('#show-history')
    .on('click', function () {
        var historyBlock =  $("#history-block");

        if (historyBlock.is(":hidden")) {
            historyBlock.show();
        } else {
            historyBlock.hide();
        }
    });


function is_form_data_changed() {
    function compare_fields(a, b) {
        if (a == b) {
            return true;
        } else if (a == null && b === "") {
            return true;
        } else if (b == null && a === "") {
            return true;
        }
        return false;
    }
    var allMatch = true;
    d3.selectAll('#form-div > :not(label)').each(function () {
        if (!compare_fields(annotation_in_db[this.name], this.value)) {
            allMatch = false;
        }
    });
    return !allMatch;
}


d3.selectAll('#form-div > :not(label)')
    .on('input', function () {
        if (is_form_data_changed()) {
            d3.selectAll('#submit-form, #reset-form').attr('disabled', null);
        } else {
            d3.selectAll('#submit-form, #reset-form').attr('disabled', 'disabled');
        }
    });


d3.select('#reset-form')
    .on('click', function () {
        d3.selectAll('#form-div > :not(label)').each(function () {
                this.value = annotation_in_db[this.name];
            });
        d3.selectAll('.form-buttons > button').attr('disabled', 'disabled');
    });


d3.select('#run-search-button')
    .on('click', function () {
        var response = $.ajax({
            type: 'POST',
            url: RUN_SEARCH_URL,
            data: {
                'tool': TOOL,
                'database': DATABASE,
                'accession': ACCESSION
            },
            dataType: 'json',
            error: function (xhr, status, error) {
                console.log(xhr.responseText);
            }
        });

        set_run_status(1);
    });


function get_selected_databases(tool) {
    var databases = [];
    $(`#${tool}-dbs.db-selection-row input`).each(
        function (i, v) {
            if ($(v).is(':checked')) {
                databases.push(v.id.slice(0,-9));
            }
        }
    );
    return databases;
}


d3.select('#submit-genome')
    .on('click', function () {
        var d = {
            'genome': NAV_ARG,
            'rerun': Boolean($('#rerun').is(':checked')),
            'tools_and_databases': {
                'blastp': get_selected_databases('blastp'),
                'rpsblast': get_selected_databases('rpsblast'),
                'hhsearch': get_selected_databases('hhsearch')
            }
        };
        var response = $.ajax({
            type: 'POST',
            url: RUN_SEARCH_FOR_PHAGE_URL,
            data: { 'data': JSON.stringify(d) },
            headers: {'X-CSRFToken': getCookie('csrftoken')},
            dataType: 'json',
            error: function (xhr, status, error) {
                console.log(xhr.responseText);
            }
        });
        location.reload();
    });


$('body').keypress(function (e) {
    if (e.key === "h" && e.ctrlKey) {
        $("#id_annotation").val('hypothetical protein');
        $("#id_flag").val(0);
        if (is_form_data_changed()) {
            d3.selectAll('#submit-form, #reset-form').attr('disabled', null);
        } else {
            d3.selectAll('#submit-form, #reset-form').attr('disabled', 'disabled');
        }
    }
});


function set_run_status(status) {

    if (TOOL === 'blastp') {
        blastp_alignment_data[DATABASE].status = status

    } else if (TOOL === 'hhsearch') {
        hhsearch_alignment_data[DATABASE].status = status

    } else if (TOOL === 'rpsblast') {
        rpsblast_alignment_data[DATABASE].status = status
    }

    if (status === 0 || status === null) {
        $('#search_glyph')
            .removeClass('glyphicon-remove-circle glyphicon-refresh')
            .addClass('glyphicon-play-circle');

        $('#run-search-button')
            .attr('disabled', (WORKER_ACTIVE ? null : 'disabled'))
            .removeClass('btn-danger running error')
            .addClass('btn-primary complete');

        $('#run-button-txt')
            .text('Run');

        if (status === 0 ) {
            var run_date = $('#search-run-date').text();
            $('#search-run-date').text(`Search ran on: ${run_date}`);

        } else {
            $('#search-run-date').text('Search never ran')
        }

    } else if (status === 1) {
        $('#search_glyph')
            .removeClass('glyphicon-remove-circle glyphicon-play-circle')
            .addClass('glyphicon-refresh');

        $('#run-search-button')
            .attr('disabled', 'disabled')
            .removeClass('btn-danger error complete')
            .addClass('btn-primary running');

        $('#run-button-txt')
            .text('Running...');

        $('#search-run-date')
            .text('Search currently in progress')

    } else if (status === 2) {
        $('#search_glyph')
            .removeClass('glyphicon-play-circle glyphicon-refresh')
            .addClass('glyphicon-remove-circle');

        $('#run-search-button')
            .attr('disabled', (WORKER_ACTIVE ? null : 'disabled'))
            .removeClass('btn-primary running complete')
            .addClass('btn-danger error');

        $('#run-button-txt')
            .text('Retry');

        var run_date = $('#search-run-date').text();
        $('#search-run-date').text(`Search failed on: ${run_date}`);
    }
}


function range(start, stop, step=1) {
    return Array(Math.ceil((stop - start) / step)).fill(start).map((x,y) => x+y*step);
}


function draw_key() {

    var key = alignment.append("g")
        .attr('id', 'key');
    var key_perc_width = 0.75;
    var rectY = 26;
    var rectHeight = 15;
    var textY = 14;

    function drawKeyHhpred() {
        //Append a linearGradient element to the defs and give it a unique id
        var defs = key.append("defs");
        var linearGradient = defs.append("linearGradient")
            .attr("id", "linear-gradient");

        linearGradient.selectAll('stop')
            .data(range(0, 101, 10))
            .enter().append("stop")
            .attr("offset", function (d) {
                return `${d}%`
            })
            .attr("stop-color", function (d) {
                return colorScale(d)
            });

        key.append('text')
            .attr('x', (CHART_WIDTH / 2))
            .attr('y', textY)
            .attr('alignment-baseline', 'middle')
            .attr('text-anchor', 'middle')
            .attr('class', 'key-title')
            .text('Probability');

        key.append('rect')
            .attr('x', CHART_WIDTH * ((1 - key_perc_width) / 2))
            .attr('y', rectY)
            .attr('width', CHART_WIDTH * key_perc_width)
            .attr('height', rectHeight)
            .style('fill', 'url(#linear-gradient)');

        key.append('text')
            .attr('x', CHART_WIDTH * ((1 - key_perc_width) / 2))
            .attr('y', textY+5)
            .text('0%');

        key.append('text')
            .attr('x', (CHART_WIDTH * key_perc_width / 2) + CHART_WIDTH / 2 - 40)
            .attr('y', textY+5)
            .text('100%')
    }

    function drawKeyBlast() {

        var keyBlockScale = d3.scaleLinear()
            .domain([0, colorScale.range().length])
            .range([
                CHART_WIDTH * ((1 - key_perc_width) / 2),
                CHART_WIDTH * ((1 - key_perc_width) / 2) + CHART_WIDTH * key_perc_width
            ]);

        function make_key_block(d, i) {
            var block = d3.create('svg:g');

            function make_text() {
                return colorScale.domain().slice().reverse()[i];
            }

            var rect = block.append('rect')
                .attr('x', keyBlockScale(i))
                .attr('y', rectY)
                .attr('width', (keyBlockScale.range()[1] - keyBlockScale.range()[0]) / colorScale.range().length)
                .attr('height', rectHeight)
                .attr('class', 'key-block')
                .style('fill', d);

            add_text_to_rect(block, rect, make_text());

            return block;
        }

        var kb = key.selectAll('.key-block')
            .data(colorScale.range().slice().reverse())
            .enter()
            .append(function (d, i) {
                return make_key_block(d, i).node();
            });

        key.append('text')
            .attr('x', (CHART_WIDTH / 2))
            .attr('y', textY)
            .attr('alignment-baseline', 'middle')
            .attr('text-anchor', 'middle')
            .attr('class', 'key-title')
            .text('E-Value');
    }

    if (TOOL === 'hhsearch') {
        drawKeyHhpred();
    } else if (TOOL === 'blastp' || TOOL === 'rpsblast') {
        drawKeyBlast();
    }
}


function add_text_to_rect(parent, rect, text) {
    // Add Text
    var t = parent.append('text');
    // http://www.w3.org/TR/AERT#color-contrast
    var rgb = rect.style('fill');
    rgb = rgb.match(/^rgb\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i);
    var o = Math.round(((parseInt(rgb[1]) * 299) +
        (parseInt(rgb[2]) * 587) +
        (parseInt(rgb[3]) * 114)) / 1000);
    var fore = (o > 125) ? 'black' : 'white';
    t.attr('x', parseFloat(rect.attr('x')) + parseFloat(rect.attr('width')) / 2);
    t.attr('y', parseFloat(rect.attr('y')) + (BAR_HEIGHT / 2));
    t.attr('alignment-baseline', 'middle');
    t.attr('text-anchor', 'middle');
    t.attr('class', 'label');
    t.style('font-size', BAR_HEIGHT * 0.8);
    t.style('fill', fore);
    t.text(text);
}


function create_subject(data) {
    var subject = d3.create('svg:g');
    subject.attr('class', 'subject');

    subject.selectAll('.hsp')
        .data(data.hits)
        .enter()
        .append(function (d, i) {
            var score_metric;
            if (TOOL === 'hhsearch') {
                score_metric = d.prob;
            } else if (TOOL === 'blastp' || TOOL === 'rpsblast') {
                score_metric = d.evalue;
            }
            return create_hsp(
                data.hit_id,
                data.row+1,
                d.query_start,
                d.query_end,
                score_metric,
                d.hit_start === 1,
                d.hit_end === data.hit_seq_len
            ).node();
        });

    // Add tooltip
    subject.attr('data-toggle', 'tooltip');
    subject.attr('title', data.hit_description);

    return subject;
}


function create_hsp(accession, row, start, stop, score, left_end_aligned, right_end_aligned) {
    var hsp = d3.create('svg:g');
    hsp.attr('class', 'hsp');

    var rect = hsp.append('rect');
    var rectY = KEY_SPACING + (ROW_SPACING + BAR_HEIGHT) * (row);
    var radius = BAR_HEIGHT / 2;
    rect.attr('height', BAR_HEIGHT);
    rect.attr('y', rectY);
    rect.attr('x', s(start) + radius);

    var rect_width = s(stop) - rect.attr('x') - radius;
    if (rect_width < 0) {
        rect_width = 0;
    }

    rect.attr('width', rect_width);

    function add_circle(left) {
        // 'left' is boolean. If not left, then right is assumed
        var circle = hsp.append('circle');
        circle.attr('r', radius);
        circle.attr('cy', rectY + radius);

        if (left) {
            circle.attr('cx', s(start) + radius);
        } else {
            circle.attr('cx', parseFloat(rect.attr('x')) + parseFloat(rect.attr('width')));
        }
        circle.style('fill', colorScale(score));
    }

    function add_broken_end(left) {
        var poly = hsp.append('polygon');
        var y1, y2, y3, y4, y5;
        y1 = parseFloat(rect.attr('y'));
        y2 = y1 + BAR_HEIGHT / 4;
        y3 = y2 + BAR_HEIGHT / 4;
        y4 = y3 + BAR_HEIGHT / 4;
        y5 = y1 + BAR_HEIGHT;

        var x1, x2;
        if (left) {
            x1 = Math.ceil(parseFloat(rect.attr('x')));
            x2 = x1 - radius;
        } else {
            x1 = Math.floor(parseFloat(rect.attr('x'))) + Math.floor(parseFloat(rect.attr('width')));
            x2 = x1 + radius;
        }

        poly.attr('points', `${x1},${y1} ${x2},${y1} ${x1},${y2} ${x2},${y3} ${x1},${y4} ${x2},${y5} ${x1},${y5}`);
        poly.style('fill', colorScale(score));
    }

    if (left_end_aligned) {
        add_circle(true);
    } else {
        add_broken_end(true);
    }

    if (right_end_aligned) {
        add_circle(false);
    } else {
        add_broken_end(false);
    }
    rect.style('fill', colorScale(score));

    add_text_to_rect(hsp, rect, accession);

    // Add click listener to jump to data in table
    hsp.on('click', function () {
        var t = $(`#${TOOL}-table-${DATABASE}`).DataTable();
        var row = t.row(`#${accession}`);
        var idx = row.index();
        var row_pos = t.rows()[0].indexOf(idx);
        var page_to_display = Math.floor(row_pos/t.page.len());
        t.page(page_to_display).draw('page');

        t.$('tr.selected').removeClass('selected');
        $(row.node()).addClass('selected');
        $('tr.selected').get(0).scrollIntoView();
    });

    return hsp;
}


function draw_query() {
    var step = 0;
    if (QUERY_LEN <= 350) {
        step = 10;
    } else if (QUERY_LEN <= 700) {
        step = 20;
    } else if (QUERY_LEN <= 1750) {
        step = 50;
    } else if (QUERY_LEN <= 3000) {
        step = 100;
    } else if (QUERY_LEN <= 6000) {
        step = 200;
    } else if (QUERY_LEN <= 10000) {
        step = 500;
    } else {
        step = 1000;
    }

    // Create and add Query to alignment
    var tick_locs = range(0, QUERY_LEN, step);
    tick_locs.push(0);
    xAxis.tickValues(tick_locs);

    var query = alignment.append('g')
        .attr('id', 'query')
        .attr('transform', `translate(0, ${KEY_SPACING})`)
        .call(xAxis);

    var rect = query.append('rect')
        .attr('y', -BAR_HEIGHT)
        .attr('x', 0)
        .attr('width', CHART_WIDTH)
        .attr('height', BAR_HEIGHT)
        .attr('id', 'queryRect');

    add_text_to_rect(query, rect, ACCESSION)
}


function determine_rows(data) {

    function overlapping(a, b) {
        return a.query_start <= b.query_end && b.query_start <= a.query_end;
    }

    var data_by_rows = [];

    for (var i = 0; i < data.length; i++) {
        for (var row = 0; row < data.length; row++) {

            if (data_by_rows.length === row) {
                data_by_rows.push([data[i]]);
                data[i].row = row;
                break

            } else {
                let has_overlap = false;
                for (var j = 0; j < data_by_rows[row].length; j++) {
                    if (overlapping(data_by_rows[row][j], data[i])) {
                        has_overlap = true;
                        break
                    }
                }
                if (!has_overlap) {
                    data_by_rows[row].push(data[i]);
                    data[i].row = row;
                    break
                }
            }
        }
    }
}


function set_svg_message(txt) {
    alignment.attr('height', KEY_SPACING + ((BAR_HEIGHT + ROW_SPACING) * 3));
        alignment.append('text')
            .attr('x', CHART_WIDTH / 2)
            .attr('y', KEY_SPACING + (ROW_SPACING + BAR_HEIGHT) * 2)
            .attr('alignment-baseline', 'middle')
            .attr('text-anchor', 'middle')
            .attr('class', 'label')
            .style('font-size', 20)
            .text(txt);
}


function draw_data(data, date, status) {
    draw_query();

    if (date !== null) {
        d3.select('#search-run-date')
            .text(date)
    }

    set_run_status(status);

    if (status === 2) {
        set_svg_message('SEARCH FAILED!');

    } else if (data === null) {
        set_svg_message('SEARCH NEVER RAN!');

    } else if (data.length > 0) {
        if (TOOL === 'hhsearch') {
            data.sort(function (a, b) {
                return b.prob - a.prob;
            });
        } else {
            data.sort(function (a, b) {
                return b.top_score - a.top_score;
            });
        }
        determine_rows(data);

        // Adjust size of SVG to fit all data
        var num_rows = data.reduce(function (accum, cv) {
            if (cv.row > accum.row) {
                return cv
            } else {
                return accum
            }
        }).row;

        var sel_num_rows = $('#svg-row-selector').children('option:selected').val();
        if (sel_num_rows !== 'All') {
            sel_num_rows = parseInt(sel_num_rows);
            if (sel_num_rows > num_rows) {
                sel_num_rows = num_rows;
            }
        } else {
            sel_num_rows = num_rows;
        }


        alignment.attr('height', KEY_SPACING + ((BAR_HEIGHT + ROW_SPACING) * (sel_num_rows + 2)));
        //KEY_SPACING + ((BAR_HEIGHT + ROW_SPACING) * (num_rows + 2)));

        alignment.selectAll('.subject')
            .data(data)
            .enter()
            .append(function (d, i) {
                return create_subject(d).node()
            });

    } else if (typeof(data) === "object") {
        set_svg_message(`NO ${TOOL.toUpperCase()} HITS`);
    }
}


function init_datatable(tool) {

    var alignment_data;
    if (tool === 'blastp') {
        alignment_data = blastp_alignment_data;
    } else if (tool === 'hhsearch') {
        alignment_data = hhsearch_alignment_data;
    } else if (tool === 'rpsblast') {
        alignment_data = rpsblast_alignment_data;
    }
    var db;

    Object.keys(alignment_data).forEach(function (db) {
        var alignments = alignment_data[db].alignments;

        if (tool === 'hhsearch') {
            var order = 4;
        } else {
            var order = 3;
        }

        var t = $(`#${tool}-table-${db}`).DataTable({
            'order': [[order, 'desc']]
        });
        t.on('click', 'div.details-control', function () {
            var tr = $(this).closest('tr');
            var hit_accession = tr.attr('id')
            var row = t.row(tr);

            if (row.child.isShown()) {
                row.child.hide();
                $(this).removeClass('shown');

            } else {
                row.child(function () {
                    var c = d3.create('pre');
                    c.text(alignments[hit_accession].text);
                    return c.node();
                }).show();
                $(this).addClass('shown');
            }
        });
    });
}


function update() {

    // Remove old data
    alignment.selectAll('*').remove();

    // Show correct table
    $('.table-div').hide();
    $(`#${TOOL}-table-div-${DATABASE}`).show();

    // Initialize color scale
    if (TOOL === 'blastp' || TOOL === 'rpsblast') {
        colorScale = d3.scaleQuantile()
            .domain([10, 1, 0.1, .001, .00001, .0000001, 0])
            .range([
                d3.interpolateMagma(0),
                d3.interpolateMagma(1/6),
                d3.interpolateMagma((1/6)*2),
                d3.interpolateMagma((1/6)*3),
                d3.interpolateMagma((1/6)*4),
                d3.interpolateMagma((1/6)*5),
                d3.interpolateMagma(1)
            ].reverse());
    } else if (TOOL === 'hhsearch') {
        colorScale = d3.scaleSequential()
            .domain([0, 100])
            .interpolator(d3.interpolateMagma);
    }
    var alignments = null;
    if (TOOL === 'blastp') {
        if (DATABASE in blastp_alignment_data && !!blastp_alignment_data[DATABASE].alignments) {
            alignments = Object.values(blastp_alignment_data[DATABASE].alignments);
        }
        draw_data(
            alignments,
            blastp_alignment_data[DATABASE].date_ran,
            blastp_alignment_data[DATABASE].status
        );
    } else if (TOOL === 'hhsearch') {
        if (DATABASE in hhsearch_alignment_data && !!hhsearch_alignment_data[DATABASE].alignments) {
            alignments = Object.values(hhsearch_alignment_data[DATABASE].alignments);
        }
        draw_data(
            alignments,
            hhsearch_alignment_data[DATABASE].date_ran,
            hhsearch_alignment_data[DATABASE].status
        );
    } else if (TOOL === 'rpsblast') {
        if (DATABASE in rpsblast_alignment_data && !!rpsblast_alignment_data[DATABASE].alignments) {
            alignments = Object.values(rpsblast_alignment_data[DATABASE].alignments);
        }
        draw_data(
            alignments,
            rpsblast_alignment_data[DATABASE].date_ran,
            rpsblast_alignment_data[DATABASE].status
        );
    }

    draw_key();
    $('[data-toggle="tooltip"]').tooltip({container: 'body'});
}


$(document).ready(function () {
    d3.select(`#${TOOL}_button`).attr('disabled', 'disabled');
    update();

    init_datatable('blastp');
    init_datatable('hhsearch');
    init_datatable('rpsblast');
    $('#history-table').DataTable();

    var csrftoken = jQuery("[name=csrfmiddlewaretoken]").val();

    function csrfSafeMethod(method) {
        return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
    }

    $.ajaxSetup({
        beforeSend: function(xhr, settings) {
            if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                xhr.setRequestHeader("X-CSRFToken", csrftoken);
            }
        }
    });
});

