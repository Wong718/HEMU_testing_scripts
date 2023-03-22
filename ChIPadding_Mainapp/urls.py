from django.urls import path, include, re_path
from Mainapp import views

urlpatterns = [
    # Home Page
    path('', views.init_scr),
    path('home', views.init_scr, name='home_page'),
    # Async task status query
    path('asynctasks/test', views.async_test, name='async_test'),
    path('asynctasks/status', views.get_async_task_progress, name='async_progress_query'),
    path('asynctasks/onhold', views.render_asyc_onhold, name='async_onhold'),
    # Site manager
    path('home/site_manage', views.site_manage, name='site_manage'),

    # Module: Transcriptome-derived analysis (fully async)
    path('gene/exp', views.gene_expression_query_async, name='gene_expression'),
    path('gene/dge', views.gene_differential_analysis_async, name='gene_DE'),
    path('gene/gokegg', views.gene_gokegg_enrichment_async, name='gene_gokegg'),
    path('gene/wgcna', views.gene_wgcna_shinyappgen_async, name='gene_wgcna'),
    path('gene/sequence/', views.gene_sequence_obtain_async, name='gene_sequence'),

    # Module: Gene family analysis
    path('genefam/hmm', views.genefam_identification_hmm_async, name='genefam_hmm'),
    path('genefam/phylo', views.genefam_phylogenetic_analysis_async, name='genefam_phylo'),
    path('genefam/htmap', views.genefam_expheatmapgen_async, name='genefam_expheatmap'),



    # Static file and gene sequence curator (async-conversion pending)
    re_path(r'^gene/static/(?P<identifier_name>\w+)/(?P<file_name>\w+.\w+);ht=(?P<frame_height>\d+);wid=('
            r'?P<frame_width>\d+)$',
            views.load_DE_staticfile),

    # Module: TE (async-conversion pending)
    path('te/exp', views.te_exp_init, name='te_expression'),
    path('te/insertion', views.te_insertion, name='te_insertion'),


    path('auxdata/userguide', views.user_guide, name='user_guide'),

    # Module: Local External Servers
    path('databrowse', views.jbrowse_catalog_render, name='jbrowse'),


    path('blast', views.seqserver, name='sequenceserver'),

]
