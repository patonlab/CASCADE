from django.urls import path, re_path
from django.conf.urls import url
from . import views
from django.views.generic import TemplateView

urlpatterns = [
    path('', views.home, name='NMR-home'),
    path('about/', views.about, name='NMR-about'),
    path('predict/', views.predict, name='NMR-predict'),
    path('sketcher/', TemplateView.as_view(template_name="cascade/sketcher.html"),
                   name='sketcher'),
    path('predict_NMR/', views.predict_NMR, name='predict_NMR'),
    path('check_task/', views.check_task, name='check_task'),
    url(r'^download/(?P<taskId>\w+)/$', views.download, name="download"),
]
