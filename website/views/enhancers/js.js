function activeAssembly(){
    var assembly = "hg19";
    if($("#assembly_mm10").is(":checked")){
        assembly = "mm10";
    }
    return assembly;
}

function updateTissues(){
    $(".selectTissues").hide();

    var assembly = activeAssembly()

    var assays = "H3K27ac";
    if($("#assaysBoth").is(":checked")){
        assays = "BothDNaseAndH3K27ac";
    } else if($("#assaysDNase").is(":checked")){
        assays = "DNase";
    }

    $("#content" + assembly + assays).show();

    // http://stackoverflow.com/a/8579673
    var aTag = $("a[href='#search']");
    $('html,body').animate({scrollTop: aTag.offset().top},
                           1400,
                           "easeOutQuint");
}

function updateExamples(){
    $(".examples").hide();
    var assembly = activeAssembly()
    $("#examples" + assembly).show();
    $("#lociSearchBox").val(annotationDefaults[assembly]["pos"]);
}

function selectAllNone(setChecked){
    // filter performance note by https://api.jquery.com/hidden-selector/
    $(".selectTissues").filter(":visible")
        .find(":checkbox").prop('checked', setChecked);
}

function select(devPt, checked){
    $(".selectTissues").filter(":visible")
        .find("*[data-selector-name='" + devPt + "']")
        .prop('checked', checked);
}

function processFormSubmitRet(event, local_url, override_data){
    event.preventDefault();

    $("#errBox").hide()

    var formData = JSON.stringify(override_data) || $("#searchForm").serializeJSON();

    //console.log(formData);

    $.ajax({
        type: "POST",
        url: local_url,
        data: formData,
        dataType: "json",
        contentType : "application/json",
        async: false, // http://stackoverflow.com/a/20235765
        success: function(got){
            if("err" in got){
                $("#errMsg").text(got["err"]);
                $("#errBox").show()
                return true;
            }

            if("url" in got){
                return window.open(got["url"], '_blank');
            }
            var w = window.open();
            $(w.document.body).html(got["html"]);
        }
    });
}

function selectIntersect(){
    $("#errBox").hide()
    $('#wait').show();

    var formData = $("#searchForm").serializeJSON();

    $.ajax({
        type: "POST",
        url: "bedsInRange",
        data: formData,
        dataType: "json",
        contentType : "application/json",
        success: function(got){
            if("err" in got){
                $("#errMsg").text(got["err"]);
                $("#errBox").show()
                return true;
            }

            if("ret" in got){
                if(null == got["ret"]){
                    $("#errMsg").text("not found");
                    $("#errBox").show()
                    return true;
                }
            }

            var assembly = activeAssembly()
            $.each(["BothDNaseAndH3K27ac", "H3K27ac", "DNase"], function(idx, assays) {
                var section = $("#content" + assembly + assays);
                section.find(":checkbox").prop('checked', false);

                $.each(got["ret"][assays], function(idx, web_id){
                    section.find(":checkbox")
                        .filter(function(){
                            return this.value == web_id})
                        .prop('checked', true);
                });
            });
            if(0 == got["ret"]["total"]){
                $("#errMsg").text("No intersecting peaks found");
                $("#errBox").show()
            }
        }
    }).always(function(){
        $('#wait').fadeOut( "slow", function() {
            // Animation complete.
        });
    });
}

function fixGallery(){
}

$(document).ready(function(){
    fixGallery();

    $(window).resize(function(){
        fixGallery();
    });

    $('[data-toggle="tooltip"]').tooltip()

    $("#selectAssembly :input").change(function() {
        updateTissues();
        updateExamples();
    });

    // http://stackoverflow.com/a/28603203
    $("#selectAssays :input").change(function() {
        updateTissues();
    });

    $("#selectAll").on('click', function(e){ e.preventDefault(); selectAllNone(true); });
    $("#selectNone").on('click', function(e){ e.preventDefault(); selectAllNone(false); });

    $.each(selectorNames, function(idx, s){
        var klass = ".select" + s;
        $(klass).on('change', function() {
            var val = this.checked;
            select(s, val);
        });
    });

    $("#hg19SNP").click(function(){
        $("#lociSearchBox").val(annotationDefaults["hg19"]["snp"]);
    });
    $("#mm10SNP").click(function(){
        $("#lociSearchBox").val(annotationDefaults["mm10"]["snp"]);
    });

    $('.tablesorter').each(function(i, obj) {
        new Tablesort(obj);
    });

    $("#searchForm").submit(function(event) {
        processFormSubmitRet(event, "ucsc");
    });

    $("#washuSubmit").click(function(event) {
        processFormSubmitRet(event, "washu");
    });

    $("#selectIntersect").click(function() {
        selectIntersect();
    });

    $("#firstGalleryPane").click(function(event){
        var fd = {"assembly":"mm10",
                  "loci":"chr19:26844000-26848000",
                  "assays":"BothDNaseAndH3K27ac",
                  "mm10BothDNaseAndH3K27ac":["neural_tube_embryonic_11_5_day"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

    $("#secondGalleryPane").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr17:40767000-40775000",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac": ["primary_natural_killer_cells_from_peripheral_blood_select",
                                              "primary_t_cells_from_peripheral_blood_select",
                                              "fetal_thymus_select"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

    $("#thirdGalleryPane").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr6:42370000-42380000",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac": ["primary_natural_killer_cells_from_peripheral_blood_select",
                                              "primary_t_cells_from_peripheral_blood_select",
                                              "fetal_thymus_select"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

    $("#fourthGalleryPane").click(function(event){
        var fd = {"assembly":"mm10",
                  "loci":"chr6:72214000-72220000",
                  "assays":"BothDNaseAndH3K27ac",
                  "mm10BothDNaseAndH3K27ac":["neural_tube_embryonic_11_5_day"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

})
