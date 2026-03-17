/*! \file   janus_opensim_audiobridge.c
 * \author OpenSim community (skeleton)
 * \copyright GNU General Public License v3
 * \brief  Janus OpenSim AudioBridge plugin (skeleton)
 *
 * This is a scaffold plugin intended as a starting point for implementing
 * listener-relative spatial mixing based on Firestorm SLData payloads (sp/sh/lp/lh).
 */

#include "plugin.h"

#include <jansson.h>

#include "../debug.h"
#include "../config.h"
#include "../utils.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Plugin information */
#define JANUS_OPENSIM_AUDIOBRIDGE_VERSION         1
#define JANUS_OPENSIM_AUDIOBRIDGE_VERSION_STRING  "0.0.1"
#define JANUS_OPENSIM_AUDIOBRIDGE_DESCRIPTION     "OpenSim-specific spatial audiobridge skeleton plugin"
#define JANUS_OPENSIM_AUDIOBRIDGE_NAME            "JANUS OpenSim AudioBridge plugin"
#define JANUS_OPENSIM_AUDIOBRIDGE_AUTHOR          "OpenSim community"
#define JANUS_OPENSIM_AUDIOBRIDGE_PACKAGE         "janus.plugin.opensim.audiobridge"

/* Plugin methods */
janus_plugin *create(void);
int janus_opensim_audiobridge_init(janus_callbacks *callback, const char *config_path);
void janus_opensim_audiobridge_destroy(void);
int janus_opensim_audiobridge_get_api_compatibility(void);
int janus_opensim_audiobridge_get_version(void);
const char *janus_opensim_audiobridge_get_version_string(void);
const char *janus_opensim_audiobridge_get_description(void);
const char *janus_opensim_audiobridge_get_name(void);
const char *janus_opensim_audiobridge_get_author(void);
const char *janus_opensim_audiobridge_get_package(void);
void janus_opensim_audiobridge_create_session(janus_plugin_session *handle, int *error);
struct janus_plugin_result *janus_opensim_audiobridge_handle_message(janus_plugin_session *handle, char *transaction, json_t *message, json_t *jsep);
json_t *janus_opensim_audiobridge_handle_admin_message(json_t *message);
void janus_opensim_audiobridge_setup_media(janus_plugin_session *handle);
void janus_opensim_audiobridge_incoming_data(janus_plugin_session *handle, janus_plugin_data *packet);
void janus_opensim_audiobridge_data_ready(janus_plugin_session *handle);
void janus_opensim_audiobridge_hangup_media(janus_plugin_session *handle);
void janus_opensim_audiobridge_destroy_session(janus_plugin_session *handle, int *error);
json_t *janus_opensim_audiobridge_query_session(janus_plugin_session *handle);

/* Plugin setup */
static janus_plugin janus_opensim_audiobridge_plugin =
    JANUS_PLUGIN_INIT(
        .init = janus_opensim_audiobridge_init,
        .destroy = janus_opensim_audiobridge_destroy,

        .get_api_compatibility = janus_opensim_audiobridge_get_api_compatibility,
        .get_version = janus_opensim_audiobridge_get_version,
        .get_version_string = janus_opensim_audiobridge_get_version_string,
        .get_description = janus_opensim_audiobridge_get_description,
        .get_name = janus_opensim_audiobridge_get_name,
        .get_author = janus_opensim_audiobridge_get_author,
        .get_package = janus_opensim_audiobridge_get_package,

        .create_session = janus_opensim_audiobridge_create_session,
        .handle_message = janus_opensim_audiobridge_handle_message,
        .handle_admin_message = janus_opensim_audiobridge_handle_admin_message,
        .setup_media = janus_opensim_audiobridge_setup_media,
        .incoming_data = janus_opensim_audiobridge_incoming_data,
        .data_ready = janus_opensim_audiobridge_data_ready,
        .hangup_media = janus_opensim_audiobridge_hangup_media,
        .destroy_session = janus_opensim_audiobridge_destroy_session,
        .query_session = janus_opensim_audiobridge_query_session,
    );

/* Plugin creator */
janus_plugin *create(void) {
    JANUS_LOG(LOG_VERB, "%s created!\n", JANUS_OPENSIM_AUDIOBRIDGE_NAME);
    return &janus_opensim_audiobridge_plugin;
}

typedef struct janus_osab_vec3 {
    double x;
    double y;
    double z;
} janus_osab_vec3;

typedef struct janus_osab_quat {
    double x;
    double y;
    double z;
    double w;
} janus_osab_quat;

typedef struct janus_opensim_audiobridge_session {
    janus_plugin_session *handle;
    gboolean datachannel_ready;
    gboolean media_ready;

    gboolean has_sp;
    gboolean has_sh;
    gboolean has_lp;
    gboolean has_lh;

    janus_osab_vec3 sp;
    janus_osab_quat sh;
    janus_osab_vec3 lp;
    janus_osab_quat lh;

    gboolean has_spatial_metrics;
    double azimuth_deg;
    int pan_0_100;
    double distance_cm;
    double gain_0_1;

    guint64 last_update_msec;
    janus_refcount ref;
} janus_opensim_audiobridge_session;

static volatile gint initialized = 0, stopping = 0;
static janus_callbacks *gateway = NULL;

static janus_mutex sessions_mutex;
static GHashTable *sessions = NULL;

static void janus_opensim_audiobridge_session_free(const janus_refcount *session_ref) {
    janus_opensim_audiobridge_session *session = janus_refcount_containerof(session_ref, janus_opensim_audiobridge_session, ref);
    g_free(session);
}

static janus_opensim_audiobridge_session *janus_opensim_audiobridge_lookup_session(janus_plugin_session *handle) {
    if(handle == NULL)
        return NULL;
    janus_mutex_lock(&sessions_mutex);
    janus_opensim_audiobridge_session *session = g_hash_table_lookup(sessions, handle);
    if(session)
        janus_refcount_increase(&session->ref);
    janus_mutex_unlock(&sessions_mutex);
    return session;
}

static double janus_osab_clamp(double v, double min_v, double max_v) {
    if(v < min_v)
        return min_v;
    if(v > max_v)
        return max_v;
    return v;
}

static double janus_osab_normalize_angle(double a) {
    while(a > M_PI)
        a -= (2.0 * M_PI);
    while(a < -M_PI)
        a += (2.0 * M_PI);
    return a;
}

static double janus_osab_quat_to_yaw_rad(const janus_osab_quat *q_in) {
    janus_osab_quat q = *q_in;

    /* Firestorm sends quaternion components multiplied by 100 over SLData. */
    if(fabs(q.x) > 2.0 || fabs(q.y) > 2.0 || fabs(q.z) > 2.0 || fabs(q.w) > 2.0) {
        q.x /= 100.0;
        q.y /= 100.0;
        q.z /= 100.0;
        q.w /= 100.0;
    }

    /* Normalize quaternion to avoid drift. */
    double qn = sqrt((q.x*q.x) + (q.y*q.y) + (q.z*q.z) + (q.w*q.w));
    if(qn > 0.0) {
        q.x /= qn;
        q.y /= qn;
        q.z /= qn;
        q.w /= qn;
    }

    /* Z-up yaw extraction. */
    return atan2(2.0 * ((q.w * q.z) + (q.x * q.y)), 1.0 - (2.0 * ((q.y * q.y) + (q.z * q.z))));
}

static void janus_osab_update_spatial_metrics(janus_opensim_audiobridge_session *session) {
    if(session == NULL)
        return;

    if(!session->has_sp || !session->has_lp || !session->has_lh) {
        session->has_spatial_metrics = FALSE;
        return;
    }

    /* World coordinates come in centimeters from Firestorm (sp/lp). */
    const double dx = session->sp.x - session->lp.x;
    const double dy = session->sp.y - session->lp.y;
    const double dz = session->sp.z - session->lp.z;

    const double distance = sqrt((dx*dx) + (dy*dy) + (dz*dz));
    const double listener_yaw = janus_osab_quat_to_yaw_rad(&session->lh);
    const double bearing = atan2(dy, dx);
    const double azimuth = janus_osab_normalize_angle(bearing - listener_yaw);

    /* -90deg(left) -> 0, 0deg(front) -> 50, +90deg(right) -> 100 */
    const double pan_f = 50.0 + (50.0 * sin(azimuth));
    const int pan = (int)llround(janus_osab_clamp(pan_f, 0.0, 100.0));

    /* Simple linear rolloff to 0 at 50m (5000cm). */
    const double gain = janus_osab_clamp(1.0 - (distance / 5000.0), 0.0, 1.0);

    session->azimuth_deg = azimuth * (180.0 / M_PI);
    session->pan_0_100 = pan;
    session->distance_cm = distance;
    session->gain_0_1 = gain;
    session->has_spatial_metrics = TRUE;
}

static gboolean janus_osab_parse_vec3(json_t *obj, janus_osab_vec3 *out) {
    if(!json_is_object(obj) || out == NULL)
        return FALSE;
    json_t *jx = json_object_get(obj, "x");
    json_t *jy = json_object_get(obj, "y");
    json_t *jz = json_object_get(obj, "z");
    if(!json_is_number(jx) || !json_is_number(jy) || !json_is_number(jz))
        return FALSE;
    out->x = json_number_value(jx);
    out->y = json_number_value(jy);
    out->z = json_number_value(jz);
    return TRUE;
}

static gboolean janus_osab_parse_quat(json_t *obj, janus_osab_quat *out) {
    if(!json_is_object(obj) || out == NULL)
        return FALSE;
    json_t *jx = json_object_get(obj, "x");
    json_t *jy = json_object_get(obj, "y");
    json_t *jz = json_object_get(obj, "z");
    json_t *jw = json_object_get(obj, "w");
    if(!json_is_number(jx) || !json_is_number(jy) || !json_is_number(jz) || !json_is_number(jw))
        return FALSE;
    out->x = json_number_value(jx);
    out->y = json_number_value(jy);
    out->z = json_number_value(jz);
    out->w = json_number_value(jw);
    return TRUE;
}

int janus_opensim_audiobridge_init(janus_callbacks *callback, const char *config_path) {
    if(g_atomic_int_get(&stopping))
        return -1;
    gateway = callback;

    janus_mutex_init(&sessions_mutex);
    sessions = g_hash_table_new(NULL, NULL);

    JANUS_LOG(LOG_INFO, "%s initialized (%s)\n", JANUS_OPENSIM_AUDIOBRIDGE_NAME,
        config_path ? config_path : "no config path");
    g_atomic_int_set(&initialized, 1);
    return 0;
}

void janus_opensim_audiobridge_destroy(void) {
    if(!g_atomic_int_get(&initialized))
        return;
    g_atomic_int_set(&stopping, 1);

    janus_mutex_lock(&sessions_mutex);
    if(sessions) {
        GHashTableIter iter;
        gpointer key = NULL;
        gpointer value = NULL;
        g_hash_table_iter_init(&iter, sessions);
        while(g_hash_table_iter_next(&iter, &key, &value)) {
            janus_opensim_audiobridge_session *session = value;
            janus_refcount_decrease(&session->ref);
        }
        g_hash_table_destroy(sessions);
        sessions = NULL;
    }
    janus_mutex_unlock(&sessions_mutex);

    janus_mutex_destroy(&sessions_mutex);

    gateway = NULL;
    g_atomic_int_set(&initialized, 0);
    g_atomic_int_set(&stopping, 0);
    JANUS_LOG(LOG_INFO, "%s destroyed\n", JANUS_OPENSIM_AUDIOBRIDGE_NAME);
}

int janus_opensim_audiobridge_get_api_compatibility(void) {
    return JANUS_PLUGIN_API_VERSION;
}

int janus_opensim_audiobridge_get_version(void) {
    return JANUS_OPENSIM_AUDIOBRIDGE_VERSION;
}

const char *janus_opensim_audiobridge_get_version_string(void) {
    return JANUS_OPENSIM_AUDIOBRIDGE_VERSION_STRING;
}

const char *janus_opensim_audiobridge_get_description(void) {
    return JANUS_OPENSIM_AUDIOBRIDGE_DESCRIPTION;
}

const char *janus_opensim_audiobridge_get_name(void) {
    return JANUS_OPENSIM_AUDIOBRIDGE_NAME;
}

const char *janus_opensim_audiobridge_get_author(void) {
    return JANUS_OPENSIM_AUDIOBRIDGE_AUTHOR;
}

const char *janus_opensim_audiobridge_get_package(void) {
    return JANUS_OPENSIM_AUDIOBRIDGE_PACKAGE;
}

void janus_opensim_audiobridge_create_session(janus_plugin_session *handle, int *error) {
    if(g_atomic_int_get(&stopping) || !g_atomic_int_get(&initialized)) {
        if(error)
            *error = -1;
        return;
    }
    if(handle == NULL) {
        if(error)
            *error = -1;
        return;
    }

    janus_opensim_audiobridge_session *session = g_malloc0(sizeof(janus_opensim_audiobridge_session));
    session->handle = handle;
    janus_refcount_init(&session->ref, janus_opensim_audiobridge_session_free);
    handle->plugin_handle = session;

    janus_mutex_lock(&sessions_mutex);
    g_hash_table_insert(sessions, handle, session);
    janus_mutex_unlock(&sessions_mutex);

    if(error)
        *error = 0;
}

struct janus_plugin_result *janus_opensim_audiobridge_handle_message(janus_plugin_session *handle, char *transaction, json_t *message, json_t *jsep) {
    (void)jsep;

    if(g_atomic_int_get(&stopping) || !g_atomic_int_get(&initialized)) {
        if(transaction)
            g_free(transaction);
        if(message)
            json_decref(message);
        return janus_plugin_result_new(JANUS_PLUGIN_ERROR, "Plugin not available", NULL);
    }

    json_t *reply = json_object();
    json_object_set_new(reply, "opensim_audiobridge", json_string("event"));

    const char *request = NULL;
    if(json_is_object(message)) {
        json_t *req = json_object_get(message, "request");
        if(json_is_string(req))
            request = json_string_value(req);
    }

    if(request && (!strcmp(request, "join") || !strcmp(request, "configure") || !strcmp(request, "leave") || !strcmp(request, "changeroom"))) {
        json_object_set_new(reply, "result", json_string("ok"));
    } else {
        json_object_set_new(reply, "error", json_string("Unsupported request in skeleton plugin"));
        json_object_set_new(reply, "error_code", json_integer(490));
    }

    if(transaction)
        g_free(transaction);
    if(message)
        json_decref(message);
    return janus_plugin_result_new(JANUS_PLUGIN_OK, NULL, reply);
}

json_t *janus_opensim_audiobridge_handle_admin_message(json_t *message) {
    (void)message;
    json_t *response = json_object();
    json_object_set_new(response, "plugin", json_string(JANUS_OPENSIM_AUDIOBRIDGE_PACKAGE));
    json_object_set_new(response, "status", json_string("ok"));
    json_object_set_new(response, "note", json_string("skeleton plugin; admin features not implemented"));
    return response;
}

void janus_opensim_audiobridge_setup_media(janus_plugin_session *handle) {
    janus_opensim_audiobridge_session *session = janus_opensim_audiobridge_lookup_session(handle);
    if(!session)
        return;
    session->media_ready = TRUE;
    janus_refcount_decrease(&session->ref);
}

void janus_opensim_audiobridge_incoming_data(janus_plugin_session *handle, janus_plugin_data *packet) {
    if(handle == NULL || packet == NULL || packet->buffer == NULL || packet->length == 0)
        return;

    janus_opensim_audiobridge_session *session = janus_opensim_audiobridge_lookup_session(handle);
    if(!session)
        return;

    char *text = g_strndup(packet->buffer, packet->length);
    if(text == NULL) {
        janus_refcount_decrease(&session->ref);
        return;
    }

    json_error_t jerr;
    json_t *root = json_loads(text, 0, &jerr);
    g_free(text);
    if(!json_is_object(root)) {
        if(root)
            json_decref(root);
        janus_refcount_decrease(&session->ref);
        return;
    }

    json_t *sp = json_object_get(root, "sp");
    json_t *sh = json_object_get(root, "sh");
    json_t *lp = json_object_get(root, "lp");
    json_t *lh = json_object_get(root, "lh");

    if(janus_osab_parse_vec3(sp, &session->sp))
        session->has_sp = TRUE;
    if(janus_osab_parse_quat(sh, &session->sh))
        session->has_sh = TRUE;
    if(janus_osab_parse_vec3(lp, &session->lp))
        session->has_lp = TRUE;
    if(janus_osab_parse_quat(lh, &session->lh))
        session->has_lh = TRUE;

    janus_osab_update_spatial_metrics(session);

    session->last_update_msec = janus_get_monotonic_time() / 1000;

    json_decref(root);
    janus_refcount_decrease(&session->ref);
}

void janus_opensim_audiobridge_data_ready(janus_plugin_session *handle) {
    janus_opensim_audiobridge_session *session = janus_opensim_audiobridge_lookup_session(handle);
    if(!session)
        return;
    session->datachannel_ready = TRUE;
    janus_refcount_decrease(&session->ref);
}

void janus_opensim_audiobridge_hangup_media(janus_plugin_session *handle) {
    janus_opensim_audiobridge_session *session = janus_opensim_audiobridge_lookup_session(handle);
    if(!session)
        return;
    session->media_ready = FALSE;
    janus_refcount_decrease(&session->ref);
}

void janus_opensim_audiobridge_destroy_session(janus_plugin_session *handle, int *error) {
    if(handle == NULL) {
        if(error)
            *error = -1;
        return;
    }

    janus_mutex_lock(&sessions_mutex);
    janus_opensim_audiobridge_session *session = g_hash_table_lookup(sessions, handle);
    if(session)
        g_hash_table_remove(sessions, handle);
    janus_mutex_unlock(&sessions_mutex);

    if(session)
        janus_refcount_decrease(&session->ref);
    handle->plugin_handle = NULL;

    if(error)
        *error = 0;
}

json_t *janus_opensim_audiobridge_query_session(janus_plugin_session *handle) {
    janus_opensim_audiobridge_session *session = janus_opensim_audiobridge_lookup_session(handle);
    if(!session)
        return NULL;

    json_t *info = json_object();
    json_object_set_new(info, "media_ready", json_boolean(session->media_ready));
    json_object_set_new(info, "datachannel_ready", json_boolean(session->datachannel_ready));
    json_object_set_new(info, "has_sp", json_boolean(session->has_sp));
    json_object_set_new(info, "has_sh", json_boolean(session->has_sh));
    json_object_set_new(info, "has_lp", json_boolean(session->has_lp));
    json_object_set_new(info, "has_lh", json_boolean(session->has_lh));
    json_object_set_new(info, "last_update_msec", json_integer((json_int_t)session->last_update_msec));
    json_object_set_new(info, "has_spatial_metrics", json_boolean(session->has_spatial_metrics));

    if(session->has_spatial_metrics) {
        json_object_set_new(info, "azimuth_deg", json_real(session->azimuth_deg));
        json_object_set_new(info, "pan_0_100", json_integer(session->pan_0_100));
        json_object_set_new(info, "distance_cm", json_real(session->distance_cm));
        json_object_set_new(info, "gain_0_1", json_real(session->gain_0_1));
    }

    if(session->has_sp) {
        json_t *sp = json_object();
        json_object_set_new(sp, "x", json_real(session->sp.x));
        json_object_set_new(sp, "y", json_real(session->sp.y));
        json_object_set_new(sp, "z", json_real(session->sp.z));
        json_object_set_new(info, "sp", sp);
    }
    if(session->has_lp) {
        json_t *lp = json_object();
        json_object_set_new(lp, "x", json_real(session->lp.x));
        json_object_set_new(lp, "y", json_real(session->lp.y));
        json_object_set_new(lp, "z", json_real(session->lp.z));
        json_object_set_new(info, "lp", lp);
    }

    janus_refcount_decrease(&session->ref);
    return info;
}
