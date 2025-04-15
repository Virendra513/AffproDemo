"""
URL configuration for AffPro project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from AffPro import views   # Load views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = [
    path('admin/', admin.site.urls),
    path('about-us/', views.aboutUS),
    path('about-us/<int:id>', views.contactDetails),
    path('', views.homePage, name='home'),
    path('index-sidebar-pl/', views.sidebar_pl, name='process_form_pl'), 
    path('index-sidebar-pp/', views.sidebar_pp, name='process_form_pp'), 
    path('index-sidebar-dd/', views.sidebar_dd, name='process_form_dd'), 
    path('index-mol-str/', views.dropdown_mol_str, name='process_form_mol_str'), 
    path('index-prot-str/', views.dropdown_prot_str, name='process_form_prot_str'), 
    path('index-li-desc/', views.dropdown_li_desc, name='process_form_li_desc'), 
    path('index-pro-desc/', views.dropdown_pro_desc, name='process_form_pro_desc'),
    path('index-pdb-read/', views.dropdown_read_pdb, name='process_form_pdb_read'),
    #path("signup/", views.signup_view, name="signup"),   # Sign-up page
    path("login/", views.login_view, name="login"),     # Login page
    path("logout/", views.logout_view, name="logout"),  # Logout functionality
    path('validate_smiles/', views.validate_smiles, name='validate_smiles'),
    path('validate_protein/', views.validate_protein_sequence, name='validate_protein'),

]



if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)