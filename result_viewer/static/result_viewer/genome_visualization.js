
$(document).ready(function() {

    function range(start, stop, step = 1) {
        return Array(Math.ceil((stop - start) / step)).fill(start).map((x,y) => x+y*step);
    }

    var feature_data = JSON.parse(document.getElementById('feature_data').textContent);

    var svg = d3.select('.chart');

    var zoom = d3.zoom()
        .on('zoom',zoomed);

    var chart = document.getElementById('chart');
    var width = chart.getBoundingClientRect().width*0.95;
    var height = chart.getBoundingClientRect().height;

    var s = d3.scaleLinear()
        .domain([0, GENOME_LENGTH])
        .range([0, width]);

    var scale_vertical = d3.scaleLinear()
        .domain([2,0])
        .range([0,50]);

    var yAxis = d3.axisLeft()
        .scale(scale_vertical)
        .tickFormat(function(d,i){
            if (i === 0){
                return '0';
            }
            else if (i%2 === 0 ){
                return '-';
            }else{
                return '+';
            }
        });

    var genome_xAxis = d3.axisBottom()
        .scale(s);

    var genome_vis;
    var features;
    var previousb = {};
    var previoust = {};

    function bottom_strand(previous, start,end){
        // console.log(previous['end']);
        if(previous['end'] !== null){
            if(start < previous.end){
                previous['end'] = end;
                if (previous['x'] === '-25'){
                    previous['x'] = '-30';
                    return '-30';
                }
                else{
                    previous['x'] = '-25';
                    return '-25';
                }
            }
            else{
                previous['end'] = end;
                previous['x'] = '-25';
                return '-25';
            }
        }else{
            previous['x'] = '-25';
            return '-25';
        }
    }

    function top_strand(previous, start, end){
        // console.log("pt: " + previous.end);
        if(previous['end'] !== null){
            if(end < previous.end ){
                previous['end'] = start;
                if (previous['x'] === '-50'){
                    previous['x'] = '-55';
                    return '-55';
                }
                else{
                    previous['x'] = '-50';
                    return '-50';
                }
            }
            else{
                previous['end'] = start;
                previous['x'] = '-50';
                return '-50';
            }
        }else{
            previous['x'] = '-50';
            return '-50';
        }
    }

    function draw() {

        var ytick_locs = range(0,3,1);
        yAxis.tickValues(ytick_locs);

        genome_vis = svg.append('g')
            .attr('id', 'genome_vis')
            .attr('transform', `translate(50, 200)`)
            .call(genome_xAxis);

        svg.append('text')
            .attr('transform','translate('+ (width/2) + ', ' + '250)')
            .style('text-anchor','middle')
            .text('Genome');

        svg.append('g')
            .attr('class', 'yAxis')
            .attr('transform', 'translate(50, 150)')
            .call(yAxis);

        svg.append('text')
            .attr('transform','rotate(-90)')
            .attr('y', 0)
            .attr('x', 0-(height/2))
            .attr('dy','1em')
            .style('text-anchor','middle')
            .text('Strand');


        features = genome_vis.selectAll('.rect')
            .data(feature_data['features'])
            .enter()
            .append('a')
            .attr('xlink:href', function(d){
                return d.href;
            })
            .append('rect').attr('class','.rect')
            .attr('x', function(d){
                if (d.strand === '+'){
                    return s(d.start);
                }else{
                    return s(d.stop);
                }
            })
            .attr('y', function(d,i){
                if (d.strand === '+'){
                    return bottom_strand(previousb, d.start, d.stop);
                }else{
                    return top_strand(previoust, d.start, d.stop);
                }
            })
            .attr('width', function(d){
                if (d.start > d.stop){
                    return s(d.start)-s(d.stop);
                }else{
                    return s(d.stop)-s(d.start);
                }
            })
            .attr('height','5')
            .style('fill', function(d){
                if (feature_data['feature_id'] == d.id)
                    return 'fuchsia';
                else if (d.flag === 'N/A' ||  d.flag === 'UNANNOTATED'){
                    return 'grey';
                }
                else if (d.flag === 'ENDOLYSIN'){
                    return 'blue';
                }
                else if (d.flag === 'REVIEW NAME'){
                    return '#C21460';
                }
                else if (d.flag === 'tRNA'){
                    return '#008080';
                }
                else if (d.flag === 'TERMINAL REPEAT'){
                    return 'purple';
                }
                else if (d.flag === 'RED'){
                    return 'red';
                }
                else if (d.flag === 'YELLOW'){
                    return 'yellow';
                }
                else if (d.annotation == 'hypothetical protein'){
                    return 'black';
                }
                else{
                    return d.flag;
                }
            })
            .style('opacity','.90')
            .attr('data-toggle', 'tooltip')
            // .attr('data-legend',function(d){return d.type;})
            .attr('title', function(d){
                return "Accession: "+ d.accession+" Feature: "+ d.start + " " + d.stop + " " + d.strand + " " + d.type + " " + d.flag +
                    ", Annotaiton: " + d.annotation + ", Pub Note: " + d.public_note + ", Private Note: " + d.private_note;
            });

    }
    draw();

    // console.log('New Js');

    function zoomed(){
        var new_x_scale = d3.event.transform.rescaleX(s);

        genome_vis.transition()
            .duration(0)
            .call(genome_xAxis.scale(new_x_scale));

        //append features
        features
            .attr('x',function(d){
               if (d.strand === '+'){
                    return new_x_scale(d.start);
                }else{
                    return new_x_scale(d.stop);
                }
            })
            .attr('width', function(d){
                if (new_x_scale(d.start) > new_x_scale(d.stop)){
                    if (new_x_scale(d.start)-new_x_scale(d.stop) >= 0){
                        // console.log(new_x_scale(d.start-d.stop));
                        return new_x_scale(d.start)-new_x_scale(d.stop);
                    }
                }else{
                    if (new_x_scale(d.stop)-new_x_scale(d.start) >= 0){
                        // console.log(new_x_scale(d.stop-d.start));
                        return new_x_scale(d.stop)-new_x_scale(d.start);
                    }
                }
            });
    }
    $('[data-toggle="tooltip"]').tooltip({container: 'body'});
    svg.call(zoom);

    var legend = svg.append('g')
        .attr('class','legend')
        .attr('transform','translate(100,20)')
        .attr('width','1000')
        .attr('height','1000');

    legend.append('circle')
        .attr('class','protein')
        .attr('cx', 0)
        .attr('cy','-5')
        .attr('r','5')
        .style('fill','black');

    legend.append('text')
        .attr('class','protein-text')
        .attr('x', '10')
        .attr('y','0')
        .style('fill','black')
        .text('Color represented by flag; black is hypothetical protein');

    legend.append('circle')
        .attr('class','tRNA')
        .attr('cx',0)
        .attr('cy',10)
        .attr('r','5')
        .style('fill','#008080');

    legend.append('text')
        .attr('class','tRNA-text')
        .attr('x', '10')
        .attr('y','15')
        .style('fill','#008080')
        .text('tRNA');

    legend.append('circle')
        .attr('class','repeat-region')
        .attr('cx',0)
        .attr('cy',25)
        .attr('r','5')
        .style('fill','purple');

    legend.append('text')
        .attr('class','repeat-region-text')
        .attr('x', '10')
        .attr('y','30')
        .style('fill','purple')
        .text('Repeat Region');

    legend.append('circle')
        .attr('class','other')
        .attr('cx',0)
        .attr('cy',40)
        .attr('r','5')
        .style('fill','grey');

    legend.append('text')
        .attr('class','other')
        .attr('x', '10')
        .attr('y','45')
        .style('fill','grey')
        .text('Unannotated, N/A');

    legend.append('circle')
        .attr('class','other')
        .attr('cx',0)
        .attr('cy',55)
        .attr('r','5')
        .style('fill','blue');

    legend.append('text')
        .attr('class','other')
        .attr('x', '10')
        .attr('y','60')
        .style('fill','blue')
        .text('Endolysin');

    legend.append('circle')
        .attr('class','other')
        .attr('cx',0)
        .attr('cy',70)
        .attr('r','5')
        .style('fill','#C21460');

    legend.append('text')
        .attr('class','other')
        .attr('x', '10')
        .attr('y','75')
        .style('fill','#C21460')
        .text('Review Name');

    legend.append('circle')
        .attr('class','other')
        .attr('cx',0)
        .attr('cy',85)
        .attr('r','5')
        .style('fill','fuchsia');

    legend.append('text')
        .attr('class','other')
        .attr('x', '10')
        .attr('y','90')
        .style('fill','fuchsia')
        .text('Current feature');

});


